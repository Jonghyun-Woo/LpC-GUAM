# 시뮬레이션 실행 (원본 Simulink)

> 루트 `CLAUDE.md`에서 분리된 상세 문서. 원본 Simulink 시뮬레이션을 실행·검증할 때 참조한다.

CLI는 없으며, 모든 것은 저장소 루트를 작업 디렉토리로 하여 MATLAB 명령창에서 실행한다.

- `RUNME.m` — 최상위 대화형 진입점. 데모 궤적(1-5)을 고르도록 사용자에게 물은 뒤 `sim(model)`과
  `simPlots_GUAM`을 호출하여 결과를 그린다.
- `setupPath.m` — 가장 먼저 실행해야 한다(직접 또는 위 스크립트를 통해). `restoredefaultpath`를 호출하고
  `ClassDef/`, `Environment/`, `lib/`(+하위 폴더), `utilities/`, `setup/`, `Exec_Scripts/`, `vehicles/`,
  `Bez_Functions/`, `Challenge_Problems/`를 MATLAB 경로에 추가한다.
- `simSetup.m` — 핵심 셋업 파이프라인: `userStruct.variants` → `SimIn`(via `setupTypes`)을 구성하고,
  스위치(`setupSwitches`)를 적용하며, 원하는 기준 궤적을 `target.RefInput`으로 해석한 후,
  `setup(SimIn, target)`, `setupParameters`, `setupVariants`, `setupBuses`를 호출하고 사용자 정의 SimOut
  버스를 구축한다.
- 데모 궤적 스크립트는 `Exec_Scripts/`에 있다(`_traj.m`로 끝나는 파일들, 그리고 `exam_RAMP.m`,
  `exam_Bezier.m`). 각 스크립트는 `userStruct`/`target`을 설정하고 `simSetup`을 호출한 뒤, 사용자가 열린
  Simulink 모델을 실행한다.
- `Challenge_Problems/RUNME.m`은 자율비행 연구(충돌 회피, 구동기 고장 등)를 위한 자기(own-ship) 궤적 +
  고장 시나리오의 별도 데모다. 그 안의 `Generate_*.m` 스크립트는 번들 데이터셋을 재생성한다.

전체 셋업 스크립트 없이 변형(variant)만 프로그래밍적으로 바꾸려면, 예:
```matlab
userStruct.variants.actType = 3; simSetup
```
변형 열거형은 `ClassDef/*Enum.m`에 있다(예: `ActuatorEnum`, `ForceMomentEnum`, `RefInputEnum`). 기본 변형
선택은 비행체별로 `vehicles/Lift+Cruise/setup/setupDefaultChoices.m`에 있다.

이 저장소에는 자동화된 테스트나 린터가 없다 — 검증은 Simulink 모델을 실행하고 `logsout{1}.Values`(보통
`SimOut`에 할당)를 살펴보거나 `vehicles/Lift+Cruise/Utils/simPlots_GUAM.m`, `utilities/Animate_SimOut.m`을
사용하여 수행한다.
