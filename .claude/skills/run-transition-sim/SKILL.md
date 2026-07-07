---
name: run-transition-sim
description: Run and validate the Refactoring/ m-code Hover2Cruise closed-loop transition simulation. Use when asked to run, execute, smoke-test, reproduce, or verify the m-code GUAM sim (Refactoring/run_transition_sim.m or LpC_GUAM().run()).
---

# run-transition-sim

## 목적
`Refactoring/` m-code 포팅의 Hover2Cruise 천이(40 s) 폐루프 시뮬레이션을 실행하고, 결과가 정상 거동인지
판정한다. 원본 Simulink(`RUNME.m`)가 아니라 m-code 포팅 전용이다.

## 실행 방법 (headless 우선)
`matlab -batch` 헤드리스 실행을 우선한다(시작 ~30 s 소요). 사용자에게 요청하지 말고 직접 실행한다.
설치된 MATLAB 버전은 머신마다 다를 수 있으므로(R2024b 또는 R2025 계열) 버전에 의존하는 코드는 피한다.
저장소 루트에서:

    matlab -batch "cd('Refactoring'); run_transition_sim"

프로그래밍적 진입(경로가 이미 잡힌 경우):

    guam = LpC_GUAM(); out = guam.run();

`matlab` CLI가 없으면(이 프로젝트 기본 가정), 사용자에게 MATLAB 명령창에서 위를 실행하고 콘솔 출력/플롯을
공유해 달라고 요청한다. 세션에서는 `!` 프리픽스로 명령을 직접 실행하도록 안내할 수 있다.

## 성공 판정 기준
- 40 s 동안 발산(NaN/Inf, 상태 폭주) 없이 안정.
- 속도 기준(`ref_vel`)을 추종.
- 북쪽 위치에 최대 ~40 ft 오버슈트가 나타나는 것은 **알려진 정상 특성(버그 아님)** — 원본 예제의 기준
  위치(구간별 선형)와 속도(램프)가 상호 불일치하고 외루프 위치 이득이 약하기(0.1) 때문. 원본 시뮬레이션
  거동과 일치한다.

## 실패 시 점검 순서
1. **경로**: `run_transition_sim.m`이 `GUAM/ controller/ config/ utils/`를 경로에 추가하는지.
2. **트림표 로드**: `controller/trim_table_Poly_ConcatVer4p0.mat`가 로드되는지.
3. **발산**: NaN/Inf면 적분 스텝(`config/SimConfig.m`)과 공력 clamp(`AeroPolynomial/check_fac_limits.m`)
   확인. 유효 비행 영역 밖 입력을 clamp하지 못하면 폭주 가능.
4. **구동기 포화**: 명령이 한계에 붙어 추종 실패면 `GUAM/EngineDynamics.m` / `GUAM/SurfaceDynamics.m`의
   rate/position 제한 확인.

## 참고
- 아키텍처/파일 맵: `Refactoring/CLAUDE.md`
- 이론(프레임/동역학/제어): `Refactoring/docs/WIKI.md`
