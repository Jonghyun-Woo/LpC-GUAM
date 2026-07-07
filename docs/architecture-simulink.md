# 아키텍처 (원본/Simulink 시뮬레이션)

> 루트 `CLAUDE.md`에서 분리된 상세 문서. 원본 Simulink 시뮬레이션의 구조를 이해·수정할 때 참조한다.

- **세 개의 핵심 구조체**가 시뮬레이션 전체를 관통한다: `SimIn`(고정 시뮬레이션 입력 / 변형 선택),
  `SimPar`(조정 가능한 시뮬레이션 파라미터), `SimOut`(로깅 출력, `logsout{1}.Values`에서 유래).
- **변형(variant) 시스템**: `userStruct.variants`가 교체 가능한 서브시스템 구현 중에서 선택한다(비행체 타입,
  기준입력 타입, 대기(atmosphere), 난류, 제어기, 구동기, 추진, 힘/모멘트 모델, EOM, 센서, 실험 타입).
  각 축(axis)마다 `ClassDef/<X>Enum.m`과 `setup/select<X>Type.m` / `utilities/variant_definitions/setup<X>Variants.m`
  쌍이 있어 Simulink 변형 서브시스템을 연결한다.
- **두 개의 교체 가능한 공력-추진 모델**(`fmType`): `Polynomial`(기본 — 빠름, 반응표면/CFD 피팅,
  `vehicles/Lift+Cruise/AeroProp/Polynomial/`에 있으며 설계상 유효 비행 영역 밖에서는 에러를 냄)과
  `SFunction`(느림, 1차 원리 스트립 이론 클래스/객체 모델, `vehicles/Lift+Cruise/AeroProp/SFunction/` 아래에
  있으며 로터/조종면 배치 구성 가능).
- **비행체 정의**는 `vehicles/Lift+Cruise/` 아래에 있다 — `ClassDef/`(비행체별 클래스), `Control/`(baseline
  LQRi 이득 스케줄 제어기, hover/transition/cruise를 통합, heading-rotated NED 프레임에서 동작),
  `Subsystems/`(`.slx` 블록: AeroProp, GearSystem, PowerSystem, Sensors, SurfEng, VehicleEOM),
  `Trim/`(오프라인 트림 루틴; 최상위 `trim_helix.m`, 비용 `mycost.m`, 제약 `nlinCon_helix.m`),
  `setup/`(비행체별 셋업, `User_SimOut/` 사용자 출력 스크립트 포함).
- **기준 궤적**: `userStruct.variants.refInputType`로 선택한다(`FOUR_RAMP`, `ONE_RAMP`, `TIMESERIES`,
  `BEZIER`, `DEFAULT`/도블릿). Bezier 곡선 생성/평가 유틸리티는 `Bez_Functions/`에 있다(`genPWCurve.m`,
  `evalPWCurve.m`, `evalSegments.m` 등). `pwcurve.waypoints` / `pwcurve.time_wpts` `.mat` 파일 또는
  `target.RefInput.Bezier` 구조체가 NED 프레임으로 표현된 웨이포인트를 공급할 수 있다.
- **신호 버스(bus)**: `utilities/bus_definitions/*.m`에 정의되고 `setup/setupBuses.m` /
  `utilities/buildBusObject.m`가 조립하며, Earth-Centered Inertial, ECEF, NED, Navigation, Velocity, Wind,
  Stability, Body 기준 프레임을 지원한다.
- **쿼터니언/회전 수학** 유틸리티(`Qmult`, `QmultSeq`, `Qinvert`, `Qtrans`, `QrotX/Y/Z`, `qcheck4`)는
  `lib/utilities/`에 있다.
- **이득 스케줄링(gain scheduling)**: `vehicles/Lift+Cruise/Control/ctrl_scheduler_GUAM.m`이 최상위
  스크립트다. 종축(longitudinal)과 횡축(lateral)을 `get_lin_dynamics_heading.m`(비행 조건에서 선형화)과
  축별 스크립트(`get_lat_dynamics_heading.m`, `ctrl_lat.m` 등)를 통해 별도로 스케줄링한다.
