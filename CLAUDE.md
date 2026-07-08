# CLAUDE.md

이 파일은 이 저장소에서 작업할 때 Claude Code (claude.ai/code)에게 지침을 제공합니다. **항상 로드되는
전역 컨텍스트**이므로 얇게 유지하고, 상세 지식은 아래 인덱스가 가리키는 문서에서 필요할 때 로드한다.

## ⚠️ 중요 규칙 (반드시 준수)

- **Discussion 및 문서 작성(`.md`, 주석 외 설명, 커밋 메시지 본문의 설명 등)은 한국어로 작성한다.**
- **코드 및 코드 내 주석(inline/block comment) 작성은 영어로 작성한다.**
- 즉, 사람에게 설명·논의하는 산문은 한국어, 소스코드와 그 안의 주석은 영어로 통일한다.
- **비자명한 변경은 계획 우선**: 탐색(`docs/` 설계 문서 참조) → 논의 → 편집 가능한 계획 파일
  (`plans/<slug>.md`) 작성 → 사용자 검토·승인 → 실행 순서를 따른다. 이 절차는 `/plan` 커맨드로
  캡슐화되어 있다. 사소한 변경은 예외.

## 프로젝트 개요

GUAM (Generic UAM Simulation)은 NASA LaRC가 개발한 범용 Lift+Cruise eVTOL 천이(transition) 비행체의
MATLAB/Simulink 시뮬레이션이다. MATLAB 셋업 스크립트로 구동되는 Simulink 모델(`GUAM.slx`)이며, 6-DOF 강체,
교체 가능한(swappable) 공력-추진/구동기/센서/제어 변형(variant) 서브시스템, 그리고 기준 궤적 생성
(램프, 시계열, 구간별 Bezier 곡선, 도블릿)을 포함한다. 빌드 시스템·패키지 매니저·테스트 프레임워크는 없으며,
MATLAB/Simulink 안에서 실행하도록 만들어진 스크립트/모델 코드다.

이 Simulink 시뮬레이션을 `Refactoring/` 아래의 순수 객체지향(OOP) MATLAB 구현으로 포팅한다.
해당 포팅은 hover-to-cruise 천이 시나리오에 대한 **실행 가능한 end-to-end 폐루프(closed-loop)
시뮬레이션**이며, 평평한 지구(flat-earth) 형식으로 제한된다(NED/Body/Aero 프레임만 사용 — ECI/ECEF 없음,
난류(turbulence) 없음). 진입점은 `Refactoring/run_transition_sim.m`이고, 전체 이론 문서는
`Refactoring/docs/WIKI.md`에 있다.

## 문서 인덱스 (작업 유형별 라우팅)

작업 대상에 따라 아래 문서를 참조한다. 상세 내용은 항상 로드하지 않고 필요할 때 해당 파일을 읽는다.

| 작업 대상 | 참조 문서 |
|---|---|
| **계획 수립 시** 기획 의도(무엇을·왜) | [`docs/PRD.md`](docs/PRD.md) |
| **계획 수립 시** 포팅 구조/클래스 관계 | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| **계획 수립 시** 과거 설계 결정 로그 | [`docs/ADR.md`](docs/ADR.md) |
| 원본 Simulink 시뮬레이션 **실행/검증** | [`docs/running-simulink.md`](docs/running-simulink.md) |
| 원본 Simulink **아키텍처/구조 수정** | [`docs/architecture-simulink.md`](docs/architecture-simulink.md) |
| `Refactoring/` m-code 포팅 작업 | [`Refactoring/CLAUDE.md`](Refactoring/CLAUDE.md) (해당 폴더 작업 시 자동 로드) |
| m-code 이론(프레임/동역학/제어) 심화 | [`Refactoring/docs/WIKI.md`](Refactoring/docs/WIKI.md) |
| Challenge Problems (자율비행/고장 시나리오) | [`Challenge_Problems/README.md`](Challenge_Problems/README.md) |
