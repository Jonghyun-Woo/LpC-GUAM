# ARCHITECTURE — 포팅 구조 (어떻게)

> **작성 가이드**: `Refactoring/` OOP 구현의 구조를 담는다. 무엇을·왜는 PRD.md, 개별 결정의
> 근거는 ADR.md를 참조한다. 세부 이론(프레임/동역학/제어 수식)은 중복 서술하지 말고
> `Refactoring/docs/WIKI.md`를 가리킨다. `/plan` 탐색 단계에서 읽힌다.
> 아래 각 섹션의 안내 문구(`<!-- -->`)를 지우고 실제 내용으로 채운다.

## 1. 전체 개요
<!-- 진입점(Refactoring/run_transition_sim.m)에서 폐루프까지의 데이터 흐름 한눈에. -->

## 2. 디렉토리 구조
<!-- Refactoring/ 하위 폴더·핵심 파일의 역할. -->

- `trajectory/ReferenceTrajectory.m` — m-file 전용 기준 궤적 생성기(classdef, static `build(scenario, dt, T)`).
  모든 폐루프 시나리오(`climb`/`althold`)를 한 곳에 모아 dense per-step 테이블(pos/vel/chi/chidot)을 반환한다.
  Simulink·lib 무의존. (Exec_Scripts/의 Simulink 전용 궤적 생성기와 별개.)
- `config/SimConfig.m` — dt/T + `scenario` 소유. `getReferenceTrajectory()`로 위 생성기를 호출해 궤적을 로드.
- `config/RSLQRConfig.m` — 게인/한계 + 컨트롤러 이산화 스텝(dt) 전담(궤적표는 보유하지 않음).

## 3. 핵심 클래스 / 모듈과 관계
<!-- 주요 클래스(동역학/공력-추진/구동기/센서/제어/궤적)와 상호 의존 관계.
     가능하면 관계 다이어그램 또는 표로. -->

## 4. 폐루프 구성 (control loop)
<!-- 기준 궤적 → 제어기 → 구동기 → 동역학 → 센서 → 제어기로 돌아오는 루프 구성. -->

## 5. 프레임 규약
<!-- 사용하는 좌표계(NED/Body/Aero)와 변환 지점. 상세 수식은 WIKI.md 링크. -->

## 6. 확장 지점 (variant/swap)
<!-- 교체 가능한 서브시스템(공력-추진, 제어 변형 등)을 어디서 어떻게 갈아끼우는가. -->

- **기준 궤적 시나리오**: `LpC_GUAM('climb')` 또는 `LpC_GUAM('althold')`로 선택(기본 `althold`).
  진입 스크립트는 `run_transition_sim.m` 상단 `scenario` 변수로 전환. 새 시나리오는
  `trajectory/ReferenceTrajectory.m`의 `build` 내 `switch`에 한 케이스를 추가하면 된다(게인/제어 코드 무수정).
