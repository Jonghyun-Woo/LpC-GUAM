# ADR — 아키텍처 결정 로그

> **작성 가이드**: 되돌리기 어렵거나 논쟁적인 설계 결정을 시간순으로 누적 기록한다.
> 하나의 결정 = 하나의 ADR 항목. 이미 내려진 결정은 수정하지 말고, 뒤집히면 새 항목에서
> "supersedes ADR-NNNN"으로 대체 사실을 남긴다. `/plan` 탐색 단계에서 과거 결정 근거로 읽힌다.
> 아래 템플릿을 복사해 항목을 추가한다.

---

## ADR-0001: 기준 궤적 생성기 분리 + 시나리오 선택

- **상태**: 채택됨
- **날짜**: 2026-07-08
- **맥락**: CBF 제어기 설계에 앞서 크루즈 고도유지 시나리오가 필요했다. 기존 Hover2Cruise 궤적은
  크루즈 진입(t=20s, 80ft) 이후에도 고도를 100ft까지 계속 올렸고(기준 위치와 w속도 불일치), 궤적이
  `RSLQRConfig.m`에 하드코딩돼 있어 시나리오를 선택할 수 없었다. 원본 궤적은 보존해야 했고,
  `Exec_Scripts/`의 궤적 생성기들은 Simulink 전용이라 m-file 전용 폐루프 포팅에서는 직접 쓸 수 없다.
- **결정**: m-file 전용 궤적 생성기 `Refactoring/trajectory/ReferenceTrajectory.m`(classdef, static
  `build(scenario, dt, T)`)를 신설해 모든 시나리오(`climb`/`althold`)를 한 곳에 모았다. `SimConfig`가
  `scenario` 프로퍼티와 `getReferenceTrajectory()`로 궤적을 로드한다. `LpC_GUAM(scenario)` 생성자로
  선택하며 기본값은 `althold`. `RSLQRConfig`는 게인/한계 + 컨트롤러 이산화 스텝(`dt`) 전담으로 정리하고,
  궤적표(`ref_*`)와 `traj_len`을 제거했다.
- **근거 / 대안**:
  - Exec_Scripts 미러 스크립트 추가안 → 기각. Simulink 전용이라 헤드리스 자동 검증이 어렵고, 우리는
    Simulink를 전혀 쓰지 않는다.
  - `RSLQRConfig` 안에 static method로 두는 안 → 기각. 게인 담당 Config에 시나리오·궤적이 섞인다.
  - Constant 이중 테이블(climb/althold 둘 다 상수로) → 기각. 단일 소스 원칙 위배, 런타임 선택 불가.
- **결과 / 트레이드오프**: 궤적 추가는 `ReferenceTrajectory.build`의 `switch` 한 곳에서 끝난다. 게인과
  궤적의 관심사가 분리됐다. `climb`은 이전 하드코딩 표를 비트 동일 재현(검증됨). 단, 컨트롤러 적분
  스텝 `dt`는 `RSLQRConfig`(컨트롤러가 `obj.cfg.dt`로 소비)와 `SimConfig` 양쪽에 남아 값(0.01)이
  중복된다 — 컨트롤러가 `SimConfig`를 보지 않기 때문이며, 두 곳을 함께 유지해야 한다.

## ADR-0002: 종방향 Liveness 필터 (HJ Reachability)

- **상태**: 채택됨 (필터 구현·수식 검증 완료 / 폐루프 개입 검증은 보류)
- **날짜**: 2026-07-08
- **맥락**: RSLQR nominal 제어기 뒷단에서 동작하는 liveness 필터가 필요했다. helperOC로 트림 조건별
  종방향 BRT 값함수 `V(x)`(WH3 하강, UH 1–20)를 계산해 보유 중이다. 논문
  *"On Safety and Liveness Filtering Using HJ Reachability"* §IV의 필터를 포팅한다. 상세는
  [`docs/refs/liveness-filter.md`](refs/liveness-filter.md).
- **결정**:
  - **필터 2종**: Least-Restrictive(eq.15, banded 경계 `V≥−eps_band`, `eps_band=1e-3`)와
    Smooth Blending(eq.20, CBF-형 QP, `gamma=5.0`) 모두 구현. 논문 권고는 높은 γ의 smooth blending.
  - **채널 범위**: **종방향(LON)만 RSLQR에 배선**. 코드는 채널-파라미터(`spec`)로 LAT도 지원하나
    미배선(LAT BRT 미보유).
  - **사영 공간**: effector-섭동 `act_lon(1:11)`(scale_factor 리미팅 **이후**)을 필터링. 가상 자세
    `act_lon(12)=theta_v`는 제외. 삽입 위치는 `pseudo_alloc`의 `u_lon` 추출 직전(최소 침습).
  - **제어-affine 특수화**: converged V(`DtV=0`) + `lon_dMax=0`(`d*=0`) → `alpha=∇V·(A_lon x)`,
    `beta=(∇V·B_lon)'`. `A/B`는 RSLQR `ctrl_lon.Ap/Bp`.
  - **QP**: Optimization Toolbox 비의존 box+halfspace 이중 이분법(`solve_box_halfspace_qp`).
  - **입력 bound**: helperOC `@GUAM_LON/init_input_bounds`와 bit-동일(BRT 제어권한 일치).
  - **스케줄링**: `tables/BRT/`의 실존 파일에서 `(UH,WH)` 인덱스 자동 파싱, UH 보간, 커버리지 밖
    pass-through. 단위 ft/rad(helperOC 시각화의 m/deg 스케일 미적용).
  - **커버리지 한계**: LON BRT는 WH3(하강)뿐 → althold(WH1–2)에선 필터 상시 pass-through(의도된 한계).
- **근거 / 대안**:
  - 모멘트공간(`mdes`) 필터 → 기각. BRT는 effector 공간에서 정의됨.
  - `quadprog` 의존 → 기각. 프로젝트 무의존 원칙(전용 이분법으로 충분, quadprog 대조로 검증).
  - 라이브러리 CBF/DeepReach → 기각. 격자 기반 BRT 값함수만 사용.
  - 채널별 코드 복제 → 기각. `spec` 파라미터화로 LON/LAT 공용.
- **결과 / 트레이드오프**:
  - Step 1–6 구현·검증 완료: `verify_liveness_vf/bounds/qp/integration.m` 모두 PASS. althold에서
    필터 on/off가 bit-동일(활성 0회)이라 회귀 안전.
  - **폐루프 개입 검증(Step 6.5)은 보류 / 원인 규명 완료(2026-07-09)**. 규명 결과 근본 원인은
    **선형 트림 BRT가 실제 비선형 GUAM 플랜트에 대한 유효한 liveness certificate가 아니라는 점**이며,
    필터/BRT/QP 구현 자체는 정확하다(선형모델 위에선 V 유지·필터 정상 작동 확인). 규명 사슬:
    - **프레임 오기(부차)**: 스케줄 축 `WH`는 heading-frame이 아니라 **body-frame w**다(트림표에서
      `WH`=`XU0(3)`=body w 정확히 일치; RSLQR.m:16-17 주석 오기). 따라서 `WH3=+11.67`은 호버에선 수직하강,
      크루즈(θ≈5°)에선 실제 관성 하강률≈0(수평비행). 이전 `lon_brt_verify`의 `pos(3)`가 body-w를 관성하강으로
      오인 적분(프레임 비정합)했으나 검증 구간(θ 작음)에선 영향 미미.
    - **발산 트리거 = pitch-rate 섭동**: 성분 격리 결과 유일 발산원은 q. `dq=+0.05` 회복, `dq=+0.15` 발산
      (θ→23°, w→455 ft/s). 저속 WH3에서 baseline pitch-rate 흡인영역 ~0.1 rad/s.
    - **필터가 못 막는 기전**: 선형 예측 `dV/dt=α+β'u`가 실제 비선형과 **부호 반대**. `∇V`의 θ성분 가중치
      지배적이고 `dθ/dt=q`는 순수 기구학(상대차수 2)이라 제어가 직접 못 건드림 → q>0면 `dV/dt`가 양수로 강제.
      V가 깊은 내부라 CBF 여유(−γV) 커서 필터 무개입 → 실제 V 상승 → 경계 도달 시 q가 이미 커 못 잡고 탈출.
    - **모델 일치 확인**: @GUAM_LON.m:40-41 — BRT의 A_lon/B_lon은 RSLQR과 동일 트림표 `Ap/Bp_lon_interp` 출처
      (필터모델=BRT모델). **Airtight 확정**: 같은 상태를 선형모델 `ẋ=A_lon·x+B_lon·u` 위에서 blend로 돌리면
      V가 −0.0007로 유지(388/1500 스텝 개입) → 실패는 순전히 linear-vs-nonlinear 간극.
    - **해결 방향(별도 작업)**: (a) 비선형 reachability(DeepReach 등)로 실제 플랜트 정합 BRT 재생성, 또는
      (b) 필터가 런타임 실제 `∇V·f(x,u)`(비선형 야코비안)를 쓰도록 확장. 현실적 천이 검증엔 런타임 외란
      주입 장치도 필요(향후 과제).

<!-- 다음 항목은 아래 템플릿을 복사해 추가한다.
## ADR-000N: <결정 제목>
- **상태**: 제안됨 | 채택됨 | 폐기됨 | 대체됨(→ ADR-NNNN)
- **날짜**: YYYY-MM-DD
- **맥락**: <어떤 문제·제약 때문에 결정이 필요했는가.>
- **결정**: <무엇을 선택했는가.>
- **근거 / 대안**: <검토한 대안과 이걸 고른 이유.>
- **결과 / 트레이드오프**: <이 결정이 낳는 영향, 감수하는 비용.>

     남길 만한 후보 결정들: flat-earth 프레임 한정 / Simulink→OOP 포팅 방침 / 제어기 선택(RSLQR·CBF). -->
