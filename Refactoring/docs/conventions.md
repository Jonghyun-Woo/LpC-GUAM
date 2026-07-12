# 규약 치트시트 (Conventions Cheat-Sheet)

> `Refactoring/` m-code에서 **부호·순서·단위 규약**을 빠르게 조회하기 위한 한 화면 요약. 동역학/공력/제어
> 코드(`RBD.m`, `Gravity.m`, `RSLQR.m`, `AeroPolynomial/`)를 만지기 전에 확인한다. 유도·상세 서술은
> [`WIKI.md`](WIKI.md)의 해당 절을 참조하고, 이 파일은 스캔용 표만 담는다. 값은 코드에서 직접 추출했다.

## 1. 프레임 (상세: WIKI §2)

- **NED 관성 프레임(flat earth)**: `rn`(북), `re`(동), `rd`(**아래가 양수**). 고도 = `-rd`.
- **Body 프레임**: 3-2-1 Euler(`phi` roll → `theta` pitch → `psi` yaw), 원점 CG. `u/v/w` = body 속도.
- **Aero 프레임**: `alpha`(받음각), `beta`(옆미끄럼각); body-axis `u < 0.01 ft/s`면 α/β = 0 처리(로깅용).
- **Heading(제어) 프레임**: 제어기가 위치/속도 오차를 heading-rotated NED에서 다룸(WIKI §2.3, §7).

## 2. 상태 벡터 (12×1, 상세: WIKI §3.1)

```
state = [rn re rd  u v w  phi theta psi  p q r]'
          1  2  3   4 5 6   7    8    9   10 11 12
```
- `1:3` NED 위치, `4:6` body 속도, `7:9` Euler 각(rad), `10:12` body 각속도.
- 유도: `dx = [p_dot; v_dot; euler_dot; omega_dot]` (`RBD.calculate_dynamics`).

## 3. 핵심 방정식·부호 (상세: WIKI §3.2)

| 항목 | 식 | 비고 |
|---|---|---|
| Body 중력 | `g_b = g·[-sθ; sφcθ; cφcθ]` | ✅ `matlab -batch` 검증(2026-07-08): θ=10°,φ=0 → `[-5.5870, 0, 31.6853]` |
| 병진 동역학 | `v_dot = -(ω × V_b) + g_b + F/m` | 중력은 `Gravity.m`으로 내부에서 더함 |
| 회전 동역학 | `ω_dot = I \ (M - ω × I·ω)` | `mldivide` 사용 |
| 위치 kinematics | `p_dot = R_b2i · V_b` | `R_b2i` = 3-2-1 DCM(body→NED) |
| Euler kinematics | `euler_dot = T_rot · ω` | `T_rot`에 `1/cθ`, `tθ` 포함(θ=±90° 특이) |

수치 적분: forward Euler, `state += dt·dx`, `dt = 0.01 s` (`config/SimConfig.m`, WIKI §3.4).

## 4. 구동기 I/O 순서 (상세: WIKI §4.2)

**트림 입력 `U0` (13×1)** — `RSLQR.interp_xu0`:
```
U0 = [flap; aileron; elevator; rudder;  lift_rotor(1:8);  pusher_rotor(9)]
       1     2        3         4        5 6 7 8 9 10 11 12  13
```

**조종면 명령 `srf_cmd` (5×1)** — `RSLQR.srface_cmd`, flaperon 믹싱:
```
[LA; RA; LE; RE; RUD]   LA = flap - aileron,  RA = flap + aileron,  LE = RE = elevator,  RUD = rudder
```
공력 입력 `d = [LA RA LE RE RUD]`(rad, 내부에서 deg 변환) — `AeroPolynomial/LpC_aero_p_v2.m`.

**엔진 명령 `eng_cmd` (9×1)** — `RSLQR.engine_cmd` = `[lift(1:8); pusher(9)] + U0(5:13)` (rad/s).
공력 입력 `n(1..9)`은 물리 로터 속도로 아래처럼 매핑됨(`LpC_aero_p_v2.m`):
```
n(1)=N1  n(2)=N3  n(3)=N5  n(4)=N7  n(5)=N2  n(6)=N4  n(7)=N6  n(8)=N8  n(9)=N9
```

## 5. 한계값 (상세: WIKI §5)

| 대상 | 한계 | 출처 |
|---|---|---|
| 조종면 (공력 유효역 clamp) | ±30° | `AeroPolynomial/LpC_interp_p_v2.m` (`check_fac_limits`) |
| 리프트 로터 (공력 유효역) | 550–1550 RPM | 동일 |
| 조종면 서보 | ±30°, `ωn = 20 rad/s`, rate 100 | `SurfaceDynamics.m` (WIKI §5.2) |
| 엔진 서보 | `ωn = 4π rad/s`, rate 100 | `EngineDynamics.m` (WIKI §5.2) |

> `check_fac_limits`는 에러를 내지 않고 경계값으로 **clamp**한다(원본과 다른 점).

## 6. 단위 (상세: WIKI §2.5)

- **기본 단위**: slug·ft·s (English). 각도는 내부적으로 rad.
- **공력 모델**은 `.deg`와 `.knot`를 소비 — `utils/Units.m`의 `u.deg`, `u.knot`로 변환.
- 단위계 전환: `Units('ft','slug')`(기본) 또는 `Units('m','kg')`(SI). English/metric 혼용 불가.
- 표준 중력 `g ≈ 32.17405 ft/s²`.
