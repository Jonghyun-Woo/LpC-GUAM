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

<!-- 다음 항목은 아래 템플릿을 복사해 추가한다.
## ADR-000N: <결정 제목>
- **상태**: 제안됨 | 채택됨 | 폐기됨 | 대체됨(→ ADR-NNNN)
- **날짜**: YYYY-MM-DD
- **맥락**: <어떤 문제·제약 때문에 결정이 필요했는가.>
- **결정**: <무엇을 선택했는가.>
- **근거 / 대안**: <검토한 대안과 이걸 고른 이유.>
- **결과 / 트레이드오프**: <이 결정이 낳는 영향, 감수하는 비용.>

     남길 만한 후보 결정들: flat-earth 프레임 한정 / Simulink→OOP 포팅 방침 / 제어기 선택(RSLQR·CBF). -->
