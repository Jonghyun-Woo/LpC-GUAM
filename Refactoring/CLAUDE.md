# CLAUDE.md — Refactoring/ (m-code 포팅)

> 이 파일은 **중첩 메모리**다. `Refactoring/` 하위 파일을 작업할 때 루트 `CLAUDE.md`에 더해 자동으로
> 컨텍스트에 로드된다. 언어 규칙 등 전역 지침은 루트 `CLAUDE.md`를 따른다.
>
> **규약 조회**: 부호·순서·단위 치트시트는 [`docs/conventions.md`](docs/conventions.md), 이론 상세는
> [`docs/WIKI.md`](docs/WIKI.md). 동역학/공력/제어 코드를 수정하기 전에 conventions를 확인한다.

## 아키텍처 (`Refactoring/` m-code 포팅)

대응하는 Simulink 블록을 순수 클래스 기반 코드로 대체한, Simulink를 사용하지 않는 병렬 OOP MATLAB
구현이다. GUAM의 기본 변형 구성(Polynomial 공력, FirstOrder 구동기, Simple EOM, baseline 제어기)을
flat-earth 가정 하에 미러링한다. 상태: **실행 가능** — MATLAB R2025a에서 Hover2Cruise 천이(40 s, 안정,
속도 기준을 추종)에 대해 헤드리스로 검증됨.

### 폴더 구조

```
Refactoring/
├── run_transition_sim.m     # 진입 스크립트 (경로 추가 → 실행 → 플롯)
├── GUAM/                     # 플랜트(plant): 동역학·공력·구동기·환경 모델
│   ├── LpC_GUAM.m            # 최상위 오케스트레이터 (폐루프 step/run)
│   ├── RBD.m                 # 6-DOF 강체 EOM (Euler 각)
│   ├── Gravity.m             # body-frame 중력 벡터
│   ├── Environment.m         # flat-earth 대기 + 정상풍(steady wind)
│   ├── AeroFrame.m           # alpha/beta/airspeed (로깅용)
│   ├── EngineDynamics.m      # 9개 로터 속도 서보 (1차 rate/pos 제한)
│   ├── SurfaceDynamics.m     # 5개 조종면 서보 (1차 rate/pos 제한)
│   ├── LpC_model_parameters.m# SACD L+C 파라미터 (참조/미러 원본)
│   └── AeroPolynomial/       # 다항식 공력-추진 모델 (원본에서 복사)
├── controller/
│   ├── RSLQR.m               # baseline 통합 RSLQR 이득스케줄 제어기 (+ LON liveness 필터 배선)
│   ├── LivenessFilter.m      # HJ-reachability liveness 필터 (LR + smooth-blend + toolbox-free QP)
│   ├── ValueFunctionLUT.m    # BRT 값함수 로드·∇V·(UH,WH,x) 스케줄 보간
│   └── trim_table_Poly_ConcatVer4p0.mat  # 트림표 + lon/lat 이득표
├── tables/                   # (gitignore) trim/ 트림표, BRT/ 종방향 BRT 값함수 (WH3, UH1..20)
├── trajectory/
│   └── ReferenceTrajectory.m # m-file 전용 궤적 생성기 (climb/althold/lon_brt_verify, static build)
├── config/
│   ├── VehicleConfig.m       # 비행체 질량/관성/기하/프롭 배치 (slug·ft)
│   ├── SimConfig.m           # dt = 0.01 s, T = 40 s, scenario 소유 + 궤적 로드
│   ├── RSLQRConfig.m         # LQR/할당 가중치 + 컨트롤러 dt (궤적표는 보유 안 함)
│   └── LivenessConfig.m      # liveness 필터 상수 + 채널 spec (helperOC BRT 정합, ft/rad)
├── utils/
│   └── Units.m               # setUnits.m 의 클래스 래퍼 (단위 변환)
└── docs/
    └── WIKI.md               # 전체 문서 (프레임/동역학/제어 이론/참고문헌)
```

### 실행 워크플로우 (per-step 파이프라인)

`run_transition_sim.m` → `guam = LpC_GUAM(scenario); out = guam.run();`(기본 `scenario='althold'`).
`LpC_GUAM`은 모든 하위 모델 객체를 생성하고, `SimConfig.getReferenceTrajectory()`로 궤적표를 로드한 뒤 첫
기준점의 트림에서 `reset()`한다. `run()`은 기준 궤적표 전체(N = 4001 스텝)를 순회하며 매 스텝 `step(ref)`를
호출하고 상태/구동기/공력프레임을 로깅한다.

`LpC_GUAM.step(ref)` 한 스텝(GUAM.slx 기본 변형과 일치, forward Euler, Δt = 0.01 s):
1. **제어**: `controller.control(state, ref)` → 엔진 명령(9×1, rad/s), 조종면 명령(5×1, rad)
2. **구동기 동역학**: `engineDynamics.step` / `surfaceDynamics.step` (rate·position 제한 1차 서보)
3. **대기**: `environment.atmosphere(-rd)` → 밀도 ρ, 음속 a
4. **공력-추진**: `R_i2b`로 body-frame 대기상대속도(`airspeed_body`) 계산 → `x_aero = [v_air; p q r]` →
   `run_LpC_aero` → body-frame 힘 `Fb`, 모멘트 `Mb`
5. **6-DOF EOM**: `rigidBody.calculate_dynamics(state, Fb, Mb)` (내부에서 중력 추가) → `dx`
6. **적분**: `state += dt * dx`, `time += dt`

### 핵심 파일

- `Refactoring/run_transition_sim.m` — 진입 스크립트: 경로 추가, 폐루프 시뮬레이션 실행, 위치/속도/자세/구동기를
  기준 대비 플롯. 프로그래밍적 사용: `guam = LpC_GUAM(); out = guam.run();`.
- `Refactoring/GUAM/LpC_GUAM.m` — 최상위 오케스트레이터(`handle` 클래스). `reset()`은 첫 기준점의 트림
  (`interp_xu0(uh0, wh0)`)에서 위치/속도/자세/각속도와 구동기 초기값을 설정한다. `step(ref)`는 위 6단계
  파이프라인을 실행하고, `run()`은 전체 시나리오를 로깅과 함께 수행한다.
- `Refactoring/GUAM/RBD.m` — flat-earth 6-DOF Euler-각 EOM(`calculate_dynamics`), 상태
  `[rn re rd u v w phi theta psi p q r]'`(12×1). 회전 동역학은 `omega_dot = I \ (M - omega × I·omega)`,
  병진은 `v_dot = -(omega × V_b) + g_b + F/m`. 중력은 `Gravity.m`(body-frame 기울임)을 통해 내부에서 더한다.
- `Refactoring/GUAM/Gravity.m` — 자세(phi, theta)로부터 body-frame 중력 벡터 `g·[-sθ; sφcθ; cφcθ]` 계산.
- `Refactoring/GUAM/AeroPolynomial/` — `vehicles/Lift+Cruise/AeroProp/Polynomial/`의 다항식 공력-추진
  모델 복사본(진입 `run_LpC_aero.m`, 블렌딩 방법 2 = sigmoid형). 조종면 입력 순서는 `[LA RA LE RE RUD]`
  (flaperon 규약). 엔진 벡터 `n`은 내부 로터에 `n(1..9) → N1 N3 N5 N7 N2 N4 N6 N8 N9`로 매핑된다.
  `check_fac_limits.m`는 여기서 에러를 내는 대신 경계값을 유지(clamp)한다.
- `Refactoring/GUAM/EngineDynamics.m` / `SurfaceDynamics.m` — ServosLib "First Order Dynamic Rate and
  Position Limited Servo"와 일치하는 1차 rate/position 제한 서보(엔진: ωn = 4π rad/s, [0, 1600/2000 RPM];
  조종면: ωn = 20 rad/s, ±30°; 둘 다 rate limit = 100). 상태는 `handle`로 스텝 간 유지.
- `Refactoring/GUAM/Environment.m` — flat-earth US Standard Atmosphere 1976 대류권(`atmosphere(h)`,
  해수면 ρ0 = 0.0023769 slug/ft³, a0 = 1116.45 ft/s) + 정상 NED 풍(`airspeed_body`); 난류 없음.
- `Refactoring/GUAM/AeroFrame.m` — 다항식 모델 규약에 따른 α/β/airspeed 계산(로깅용). body-axis u가
  임계값(0.01 ft/s) 미만이면 α/β를 0으로 처리.
- `Refactoring/controller/RSLQR.m` — baseline 통합 RSLQR 제어기 포팅(AIAA 2021-0999). `handle` 클래스로
  서보-보상기 적분기(`ctrl_lon_i`, `ctrl_lat_i`)와 가상 자세 상태(`theta_v`, `phi_v`)를 스텝 간 유지한다.
  제어 흐름(`control`):
    1. `guidance_error` — body-frame 위치 오차 `e_pos_body`와 heading 오차 `e_chi` 계산
    2. `scheduled_params(u_cmd, w_cmd)` — (UH, WH) 격자에 대한 트림(`interp_xu0`)과 lon/lat 이득행렬의
       쌍선형(bilinear) 보간(`interp_mtrx`)
    3. `perturb_lin_ctrl` — 트림 대비 섭동 상태 `Xlon/Xlat`와 명령 `Xlon_cmd/Xlat_cmd` 구성(외루프 위치
       이득 0.1 적용)
    4. `lon_ctrl` / `lat_ctrl` — 서보-보상기 적분 + LQR로 원하는 모멘트 `mdes` 산출
    5. `pseudo_alloc` — 가중 유사역행렬 할당 `M = W⁻¹Bᵀ(BW⁻¹Bᵀ)⁻¹`에 엔진/엘리베이터/플랩 한계 스케일링을
       적용, **scale_factor 리미팅 후 `act_lon(1:11)`을 `liveness_lon.filter`로 사영**(LON liveness 필터,
       커버리지 밖·mode 'off'면 pass-through), 가상 자세 효과기 θ_v/φ_v 추출, 트림 `U0`를 더해 명령 생성.
  트림 입력 순서: `U0 = [flap; ail; ele; rud; eng1..9]`. `rotm_i2b`/`rotm_b2i`는 정적(static) 회전행렬
  유틸리티(3-2-1 Euler, NED↔body).
- `Refactoring/controller/LivenessFilter.m` / `ValueFunctionLUT.m` — HJ-reachability liveness 필터
  (AIAA liveness 논문 §IV; 상세 [`docs/refs/liveness-filter.md`](../docs/refs/liveness-filter.md)).
  `LivenessFilter`는 LR(eq.15)·smooth-blending(eq.20) + toolbox-free box+halfspace QP.
  `ValueFunctionLUT`는 `tables/BRT/`의 종방향 BRT 값함수를 로드해 `(UH,WH,x)`에서 `V`·`∇V` 보간.
  **단위 ft/rad**(helperOC 시각화의 m/deg 스케일 미적용). RSLQR 생성자가 `liveness_lon`(기본 mode
  'blend')을 소유하며 LON만 배선(LAT는 코드 지원·미배선). **한계**: 현 LON BRT는 WH3(하강)뿐 → althold
  (WH1–2)에선 상시 pass-through. 폐루프 개입 검증은 보류(ADR-0002 참조: WH3는 baseline이 섭동에 불안정).
- `Refactoring/config/VehicleConfig.m` — SACD L+C 질량/관성/기하/프롭 배치(slug·ft). `LpC_model_parameters.m`을
  `Constant` 프로퍼티로 옮기고, 생성자에서 프롭 hub 회전행렬(`Prop_R_BH`)·회전축 단위벡터(`Prop_rot_axis_e`)·
  추력 방향(`p_T_e`)·역관성(`inv_I`)을 계산한다. 아래 두 방식으로 공력 모델의 `Model` 인자로 전달된다.
- `Refactoring/trajectory/ReferenceTrajectory.m` — m-file 전용 기준 궤적 생성기. static
  `build(scenario, dt, T)`가 `climb`/`althold` 시나리오의 dense 테이블(`pos/vel/chi/chidot`, N=4001)을
  반환한다. `climb`은 이전에 `RSLQRConfig`에 하드코딩됐던 표를 비트 동일 재현. Simulink·lib 무의존.
- `Refactoring/config/SimConfig.m` — Δt = 0.01 s, T = 40 s, `scenario`(기본 `althold`) 소유.
  `getReferenceTrajectory()`가 `ReferenceTrajectory.build`를 호출해 궤적을 로드한다.
- `Refactoring/config/RSLQRConfig.m` — LQR/할당 가중치(`Qlon/Rlon/Wlon`, `Qlat/Rlat/Wlat`)와 차원 상수,
  구동기 한계, 컨트롤러 서보-보상기 적분 스텝 `dt`(RSLQR가 `obj.cfg.dt`로 소비). 기준 궤적표는 더 이상
  보유하지 않는다(→ `ReferenceTrajectory`/`SimConfig`).
- `Refactoring/utils/Units.m` — `lib/utilities/setUnits.m`의 클래스 래퍼. 공력 모델은 `.deg`와 `.knot`를
  소비한다. `Units('ft','slug')`(GUAM 기본) 또는 `Units('m','kg')`(SI) 선택 가능.
- `Refactoring/docs/WIKI.md` — 전체 문서: 수정 이력, 프레임, 동역학/공력/구동기/환경 모델, 제어 이론, 논문 참고문헌.

### 참조 시나리오 (Hover2Cruise 천이, 40 s)

`ReferenceTrajectory.build(scenario, dt, T)`가 생성, `LpC_GUAM(scenario)`로 선택(기본 `althold`):
- 공통 0-20 s: 수직 상승 0 → 80 ft (초기 상승률 8 ft/s)
- 공통 20-40 s: 전진비행 15 ft/s로 가속
- `althold`(기본): 크루즈 구간 고도 80 ft 유지(기준 w=0과 정합). 폐루프 최종 Down ≈ −79 ft.
- `climb`(원본): 크루즈 구간 100 ft까지 계속 상승(기준 위치는 상승, 기준 w=0 → 위치/속도 불일치).
  폐루프 최종 Down ≈ −89 ft.

알려진 특성: 위치 추종에 유계(bounded) 지연/오버슈트(북쪽 최대 ~40 ft)가 나타난다. 이는 기준 위치(구간별
선형)와 속도(램프)가 상호 불일치하고 외루프 위치 이득이 약하기(0.1) 때문이며, 원본 시뮬레이션의 거동과
일치하는 것으로 버그가 아니다.