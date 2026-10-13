# 논문 초안 참고문헌 요약

[논문 초안](https://claude.ai/code/artifact/ad0b7359-0326-4376-b53d-5e9b832918e2)에 인용한 29편을 한 번씩 읽어 볼 수 있게, 초안에서 맡은 역할별로 묶어 요약했습니다(2026-10-14 기준).

- 논문마다 네 줄입니다.
  - **질문:** 그 논문이 푼 문제
  - **핵심:** 방법과 결과
  - **우리 논문에서:** 초안의 어느 절에서 무엇의 근거로 쓰는지
  - **읽을 곳:** 시간이 없을 때 볼 부분
- `[ ]`를 읽은 뒤 `[x]`로 바꿔 진행 상황을 기록하세요.
- **확인 상태**
  - ✅ 웹에서 서지와 초록을 확인함
  - ◐ 내용은 널리 알려졌지만 권·쪽 번호를 기억에 기댐
  - `【확인 필요】` 투고 전에 원문 확인이 꼭 필요함
- 원문 번역이 아니라 짧게 재구성한 요약이고, 수치는 확인한 것만 넣었습니다. 예전 FC→SC 방향의 논문은 [SUMMARIES.md](SUMMARIES.md)에 있습니다.

**먼저 읽을 순서 (반나절 분량)**
1. Lin et al. (2024): 가장 가까운 선행 비교. 우리가 무엇을 더했는지가 여기서 분명해집니다.
2. Pedigo et al. (2023): 두 네트워크가 "같다"를 통계적으로 어떻게 정의하는지
3. Betzel et al. (2016): 우리 생성 모델과 energy의 출처
4. Schlegel et al. (2024)와 Witvliet et al. (2021): 같은 종 안 기준선의 근거
5. Richter & Schneidman (2025): 우리 규칙 전이 실험이 메우는 빈자리
6. Puxeddu et al. (2024): "규모가 계통보다 중요하다"는 해석의 배경

---

## A. 핵심 선행 연구 — 우리 논문과 직접 비교되는 것

- [ ] **Lin, A. et al. (2024). Network statistics of the whole-brain connectome of *Drosophila*. *Nature* 634, 153–165.** ✅
  - 질문: 성체 초파리 전뇌(FlyWire)의 네트워크 통계는 어떻고, 다른 동물과 비교하면 어떤가?
  - 핵심: 선충 자웅동체·수컷, 제브라피시 후뇌, 생쥐 V1과 같은 표에 놓았습니다. 상호성과 클러스터링은 각자의 널 대비로 보면 비슷했습니다. 리치클럽은 초파리가 뉴런의 약 30%, 선충이 11개(약 4%)였습니다. 모듈성은 다루지 않았고, 크기가 다른 네트워크의 직접 비교는 조심하라고 적었습니다.
  - 우리 논문에서: Introduction과 Related work의 출발점입니다. 우리 기여는 이 비교에 노드 수·밀도·범위 맞춤과 같은 종 안 기준선, 규칙 수준 전이를 더한 것입니다.
  - 읽을 곳: 종 간 비교 표(Table 2), 리치클럽 절

- [ ] **Pedigo, B.D., Powell, M., Bridgeford, E.W., Winding, M., Priebe, C.E. & Vogelstein, J.T. (2023). Generative network modeling reveals quantitative definitions of bilateral symmetry exhibited by a whole insect brain connectome. *eLife* 12, e83739.** ✅
  - 질문: 유충 초파리 뇌의 좌우 반구는 통계적으로 "같은" 네트워크인가?
  - 핵심: ER, 확률적 블록 모델 같은 생성 모델로 두 반구를 검정했습니다. 전체 연결 확률과 세포 유형 간 연결 확률이 모두 유의하게 달랐고, 밀도를 다시 맞추거나 약한 연결을 빼면 대칭으로 볼 수 있었습니다.
  - 우리 논문에서: 좌우 반구를 같은 종 안 반복으로 쓰는 근거이자 주의점입니다(3.5절). "같다"의 판정이 보정 방식에 달려 있다는 교훈은 Discussion의 "보존을 어떻게 말해야 하나"와 이어집니다. 우리가 쓴 유충 데이터 파일도 이 논문의 공개 저장소에서 받았습니다.
  - 읽을 곳: Introduction, 밀도 보정을 다룬 결과 절

- [ ] **Schlegel, P. et al. (2024). Whole-brain annotation and multi-connectome cell typing of *Drosophila*. *Nature* 634, 139–152.** ✅
  - 질문: hemibrain과 FlyWire, 두 개체의 커넥톰은 얼마나 같은가?
  - 핵심: 세포 유형 수와 강한 연결은 대체로 안정적이었지만, 연결 가중치는 개체 안팎에서 꽤 달랐습니다. 10 시냅스 이상이거나 표적 입력의 1% 이상인 연결은 잘 보존됩니다.
  - 우리 논문에서: 초파리 쪽 개체 간 변이의 근거입니다. 성체 시냅스 위치에 쓴 FlyWire 주석 파일(세포체 좌표)도 이 논문의 공개 자료입니다.
  - 읽을 곳: hemibrain–FlyWire 비교 결과, 연결 보존 기준

- [ ] **Witvliet, D. et al. (2021). Connectomes across development reveal principles of brain maturation. *Nature* 596, 257–261.** ✅
  - 질문: 선충 뇌는 발달하며 어떻게 바뀌고, 개체마다 얼마나 다른가?
  - 핵심: 같은 유전형 8개체의 뇌(신경환)를 L1부터 성체까지 재구성했습니다. 뉴런 배치와 형태는 같았지만 연결은 개체마다 상당히 달랐고, 일부 연결은 "가변"으로 분류됐습니다.
  - 우리 논문에서: 선충 쪽 같은 종 안 기준선(8개체)이자 규칙 전이 실험의 선충 데이터입니다(3.1, 4.3, 4.4절).
  - 읽을 곳: 개체 간 변이 결과, 연결 분류(안정·가변·발달)

- [ ] **Cook, S.J., Kalinski, C.A., Loer, C.M. et al. (2025). Comparative connectomics of two distantly related nematode species reveals patterns of nervous system evolution. *Science*. doi:10.1126/science.adx2143** ✅
  - 질문: 1억 년 이상 갈라진 두 선충(*C. elegans*, *P. pacificus*)의 커넥톰은 어디까지 같은가?
  - 핵심: 상동 뉴런 단위로 비교해 보존된 핵심 커넥톰을 정의했습니다. 세포 사멸, 위치, 돌기 투사, 연결의 변화는 뇌 전체에 고르게 퍼져 있었고, 연결이 달라도 행동이 유지되는 "circuit drift"를 보고했습니다.
  - 우리 논문에서: 종 간 비교의 가장 강한 선례지만 노드 대응이 있는 같은 문 안의 비교입니다. 우리는 노드 대응이 없는 두 문 사이를 다룹니다(Related work). Limitations 7번의 확장 후보이기도 합니다.
  - 읽을 곳: 보존된 핵심 커넥톰 정의, circuit drift

- [ ] **Richter, O. & Schneidman, E. (2025). Building the connectome of a small brain with a simple stochastic developmental generative model. *PNAS* 122(47), e2504913122.** ✅ (정정 공지 있음)
  - 질문: 단순한 확률적 발달 규칙으로 선충 커넥톰을 만들 수 있는가?
  - 핵심: 세포 유형, 출생 시점, 세포체 거리, 상호성, 가지치기로 시냅스 존재와 차수·모티프 통계를 예측했습니다. 한 개체에서 배운 규칙이 다른 개체를 예측했고, 다른 종 적용은 제안만 했습니다.
  - 우리 논문에서: 규칙 전이 실험이 메우는 빈자리입니다. 세포 유형 정보가 필요해 종 간에는 쓸 수 없다는 점은 우리가 경제적 생성 모델을 고른 이유입니다(3.6절).
  - 읽을 곳: 모델 정의, 개체 간 예측 결과, 정정 공지

- [ ] **(유충 vs 성체 초파리 비교) Brain rewiring during development: a comparative analysis of larval and adult *Drosophila melanogaster* connectomes. *Netw Neurosci* 9(4), 1299 (2025).** `【확인 필요: 저자】`
  - 질문: 유충 뇌와 성체 뇌의 네트워크 구조는 어떻게 다른가?
  - 핵심: 유충만 리치클럽을 보였고, 변태 후 버섯체 뉴런의 중심성이 줄어든다고 보고했습니다.
  - 우리 논문에서: 같은 종이라도 규모·발달 단계가 다르면 구조가 다르다는 우리 결과(유충 vs 성체)와 같은 방향입니다(Related work, Discussion).
  - 읽을 곳: 리치클럽 비교 결과

## B. 데이터 출처 — 인용 의무가 있는 것

- [ ] **Dorkenwald, S. et al. (2024). Neuronal wiring diagram of an adult brain. *Nature* 634, 124–138.** ✅
  - 질문: 성체 초파리 뇌 전체의 뉴런 단위 배선도를 만들 수 있는가?
  - 핵심: FlyWire 공동체가 성체 암컷 초파리 뇌 전체(약 14만 뉴런)를 재구성해 공개했습니다.
  - 우리 논문에서: 성체 초파리 데이터(연결 표 v783, 시냅스 표)의 출처입니다(3.1절).
  - 읽을 곳: 데이터 공개 범위, 시냅스 검출 신뢰도(cleft score) 설명

- [ ] **Scheffer, L.K. et al. (2020). A connectome and analysis of the adult *Drosophila* central brain. *eLife* 9, e57443.** ◐
  - 질문: 성체 초파리 중심뇌의 큰 부분을 뉴런·시냅스 단위로 재구성할 수 있는가?
  - 핵심: hemibrain 데이터셋을 공개했습니다. 오른쪽 반구 중심이고 일부 영역은 잘려 있습니다.
  - 우리 논문에서: 회로 수준 비교(4.1절)의 hemibrain 회로 5개 출처이자, FlyWire와 다른 개체로서의 반복입니다.
  - 읽을 곳: 데이터 범위(어느 영역이 잘렸는지)

- [ ] **Winding, M., Pedigo, B.D. et al. (2023). The connectome of an insect brain. *Science* 379, eadd9330.** ✅
  - 질문: 곤충 뇌 전체의 시냅스 수준 배선도는 어떤 구조인가?
  - 핵심: 1령 유충 초파리 뇌의 뉴런 3,016개와 시냅스 548,000개를 재구성했습니다. 반복 회로와 다층 구조를 보고했습니다.
  - 우리 논문에서: 유충 데이터의 출처입니다. 선충처럼 신경계 전체가 수천 뉴런이라 범위 맞춤의 핵심입니다(4.3절).
  - 읽을 곳: 네트워크 구조 개요 그림

- [ ] **Chen, B.L., Hall, D.H. & Chklovskii, D.B. (2006). Wiring optimization can relate neuronal structure and function. *PNAS* 103, 4723–4728.** ◐
  - 질문: 선충 뉴런의 위치는 배선 길이를 최소화하도록 놓여 있는가?
  - 핵심: 뉴런 배치를 배선 비용 최소화 모델로 상당 부분 설명했고, 그 과정에서 정리한 연결 데이터가 이후 데이터셋의 바탕이 됐습니다.
  - 우리 논문에서: 선충 전신 데이터(CoMuNeLab)가 인용을 요구하는 논문이자, 거리 항(D^η)의 생물학적 근거입니다(3.6절).
  - 읽을 곳: 배선 최적화 모델과 결과 그림

- [ ] **De Domenico, M., Porter, M.A. & Arenas, A. (2015). MuxViz: a tool for multilayer analysis and visualization of networks. *J Complex Netw* 3, 159–176.** ◐
  - 질문: 다층 네트워크를 분석·시각화하는 도구를 만들 수 있는가?
  - 핵심: 다층 네트워크 분석 도구 MuxViz를 소개했습니다.
  - 우리 논문에서: 선충 다층 커넥톰 데이터셋이 인용을 요구하는 논문입니다. 내용상 비중은 없습니다.
  - 읽을 곳: 읽지 않아도 됩니다(인용만).

- [ ] **Skuhersky, M. et al. (2022). Toward a more accurate 3D atlas of *C. elegans* neurons. *BMC Bioinformatics* 23, 195.** ✅(내용) `【확인 필요: 권·번호】`
  - 질문: 개체 간 변이를 반영한 선충 뉴런의 3D 위치 아틀라스를 만들 수 있는가?
  - 핵심: NeuroPAL 데이터로 뉴런마다 최소 7마리의 위치를 정렬해 전신 3D 아틀라스를 공개했습니다.
  - 우리 논문에서: 선충 뉴런 위치(거리 항)의 출처입니다. 세포체 위치라는 한계는 Limitations 3번에 있습니다.
  - 읽을 곳: 아틀라스 구성 방법, 위치 변동 크기

## C. 생성 모델 — 규칙 전이 실험의 근거

- [ ] **Vértes, P.E. et al. (2012). Simple models of human brain functional networks. *PNAS* 109, 5868–5873.** ◐
  - 질문: 몇 개의 단순한 규칙으로 인간 뇌 네트워크의 위상을 재현할 수 있는가?
  - 핵심: 거리 비용과 위상적 이득(클러스터링)을 결합한 두 매개변수 "경제적" 모델이 fMRI 기능 네트워크의 통계를 잘 재현했습니다.
  - 우리 논문에서: 경제적 생성 모델 계열의 출발점입니다(Related work, 3.6절).
  - 읽을 곳: 모델식과 매개변수 공간 그림

- [ ] **Betzel, R.F. et al. (2016). Generative models of the human connectome. *NeuroImage* 124, 1054–1064.** ◐
  - 질문: 어떤 연결 규칙이 인간 구조 커넥톰을 가장 잘 만들어 내는가?
  - 핵심: 세 데이터셋 380명에 13개 생성 모델을 맞췄습니다. 먼 연결에 벌점을 주고 이웃이 비슷한 영역끼리 연결을 선호하는 모델(거리 + matching index)이 가장 잘 맞았습니다. 적합도는 차수·클러스터링·매개 중심성·연결 길이 분포의 KS 최댓값(energy)입니다.
  - 우리 논문에서: 모델식의 거리·이웃 유사성 항, energy 정의, 격자 탐색 방식의 출처입니다(3.6절). 우리 데이터에서는 이웃 유사성 항이 선택되지 않았다는 점이 이 논문과 다릅니다.
  - 읽을 곳: 모델 목록 표, energy 정의, 최적 모델 결과

- [ ] **Akarca, D. et al. (2021). A generative network model of neurodevelopmental diversity in structural brain organization. *Nat Commun* 12, 4216.** ◐
  - 질문: 생성 모델의 매개변수로 개인차를 설명할 수 있는가?
  - 핵심: 아동 코호트에서 개인별로 생성 모델을 맞추고, 매개변수 차이가 인지 능력 같은 개인차와 관련됨을 보였습니다.
  - 우리 논문에서: "매개변수 = 개체·종의 배선 특성 요약"이라는 관점의 근거입니다. 우리는 이것을 종 간으로 넓혔습니다(3.6절).
  - 읽을 곳: 개인별 맞추기 방법

- [ ] **Akarca, D. et al. (2023). Homophilic wiring principles underpin neuronal network topology in vitro. *eLife* 12, e85300.** `【확인 필요: 저자】`
  - 질문: 뉴런 수준에서도 경제적 생성 규칙이 성립하는가?
  - 핵심: 배양 뉴런 네트워크에서도 이웃 유사성 규칙이 공간 규칙만 쓰는 모델보다 잘 맞았습니다.
  - 우리 논문에서: 뉴런 수준 커넥톰에 경제적 모델을 쓰는 근거입니다. 우리 결과(이웃 유사성 항 미선택)와 다른 점이라 Discussion에서 다룰 거리입니다.
  - 읽을 곳: 모델 비교 결과

- [ ] **Barabási, A.-L. & Albert, R. (1999). Emergence of scaling in random networks. *Science* 286, 509–512.** ◐
  - 질문: 실제 네트워크의 허브(멱법칙 차수 분포)는 어떻게 생기는가?
  - 핵심: 네트워크가 성장하면서 새 노드가 이미 연결이 많은 노드에 붙기 쉬우면(선호적 연결) 허브가 생긴다는 모델을 제시했습니다.
  - 우리 논문에서: 모델의 β 항(허브 형성)의 근거입니다. 우리 핵심 결과인 "성체만 β가 크다"의 해석 틀입니다(3.6, 4.4절).
  - 읽을 곳: 짧은 논문이라 전체

- [ ] **"Coming up short: Generative network models fail to accurately capture long-range connectivity." *Netw Neurosci* 9(4), 1377 (2025). doi:10.1162/NETN.a.35** `【확인 필요: 저자】`
  - 질문: 생성 모델은 장거리 연결과 전체 구조를 얼마나 정확히 재현하는가?
  - 핵심: 생성 모델은 장거리 연결의 위치를 잘 맞추지 못하고, 흔히 쓰는 평가 지표도 모델과 실제의 유사성을 제대로 재지 못한다고 보고했습니다.
  - 우리 논문에서: 모델의 한계에 대한 경고입니다. 그래서 연결 길이 분포를 적합도에 넣고, 절대 적합도 대신 같은 종 안 전이에 대한 상대 비교로 판정했습니다(Related work, Limitations 5번).
  - 읽을 곳: 평가 지표 비판 부분

## D. 비교 커넥톰믹스 배경

- [ ] **van den Heuvel, M.P., Bullmore, E.T. & Sporns, O. (2016). Comparative connectomics. *Trends Cogn Sci* 20, 345–361.** ◐
  - 질문: 선충부터 인간까지 커넥톰에서 공통으로 나타나는 원리는 무엇인가?
  - 핵심: 모듈 구조, 허브와 리치클럽, 배선 비용 절약이 종을 넘어 반복된다고 정리한 리뷰입니다.
  - 우리 논문에서: "배선 원리가 보존된다"는 주장의 출처입니다. 우리는 이 보존이 방향의 보존인지 값의 보존인지를 검증합니다(Introduction).
  - 읽을 곳: 종 간 공통 원리를 정리한 그림과 표

- [ ] **Puxeddu, M.G. et al. (2024). Relation of connectome topology to brain volume across 103 mammalian species. *PLoS Biol* 22, e3002489.** ◐
  - 질문: 포유류 103종에서 커넥톰 위상은 뇌 크기와 어떤 관계인가?
  - 핵심: 같은 측정 방식(확산 MRI)으로 보면 뇌 크기에 따라 위상이 체계적으로 바뀌고, 큰 뇌일수록 공간 제약이 강해집니다.
  - 우리 논문에서: "전이는 계통보다 규모를 따라간다"는 해석의 가장 가까운 선례입니다(Introduction, Discussion).
  - 읽을 곳: 뇌 크기와 위상 지표의 관계 그림

- [ ] **Bénichou, A., Masson, J.-B. & Vestergaard, C.L. (2024). Compression-based inference of network motif sets. *PLoS Comput Biol*. doi:10.1371/journal.pcbi.1012460** ✅
  - 질문: 널 모델을 미리 고르지 않고 의미 있는 모티프를 찾을 수 있는가?
  - 핵심: 압축 기반 방법으로 멍게, 갯지렁이, 선충, 초파리(여러 발달 단계·뇌 영역) 커넥톰의 모티프를 비교했습니다.
  - 우리 논문에서: 널 모델 선택이 결론을 바꾼다는 문제를 다른 방식으로 푼 예입니다(Related work).
  - 읽을 곳: 종 간 비교 결과 그림

- [ ] **Bullmore, E. & Sporns, O. (2012). The economy of brain network organization. *Nat Rev Neurosci* 13, 336–349.** ◐
  - 질문: 뇌 네트워크는 배선 비용과 기능적 이점 사이에서 어떻게 절충하는가?
  - 핵심: 뇌는 배선 비용을 최소화하지만 완전히 최소화하지는 않고, 비용이 드는 장거리 연결과 허브로 효율을 얻는다고 정리한 리뷰입니다.
  - 우리 논문에서: 모델의 거리 항(배선 비용)과 "경제적" 모델이라는 이름의 근거입니다(3.6절).
  - 읽을 곳: 비용–효율 절충 그림

## E. 방법론 — 정규화, 널 모델, 크기·밀도, 전이 지수

- [ ] **Maslov, S. & Sneppen, K. (2002). Specificity and stability in topology of protein networks. *Science* 296, 910–913.** ◐
  - 질문: 관측된 패턴이 차수 분포만으로 설명되는가?
  - 핵심: 각 노드의 차수를 유지한 채 연결 쌍을 맞바꾸는 재배선 널 모델을 제안했습니다.
  - 우리 논문에서: 모듈성·클러스터링·리치클럽 정규화의 차수 보존 널입니다(3.3절).
  - 읽을 곳: 재배선 방법 설명

- [ ] **Colizza, V., Flammini, A., Serrano, M.A. & Vespignani, A. (2006). Detecting rich-club ordering in complex networks. *Nat Phys* 2, 110–115.** ◐
  - 질문: 리치클럽을 어떻게 올바르게 측정하는가?
  - 핵심: 허브는 우연으로도 서로 연결되기 쉬워서, 원시 리치클럽 계수는 항상 높게 나옵니다. 차수 보존 널 대비 비율로 봐야 한다고 보였습니다.
  - 우리 논문에서: 리치클럽 정규화의 근거입니다(3.3절).
  - 읽을 곳: 정규화 정의와 예시

- [ ] **van Wijk, B.C.M., Stam, C.J. & Daffertshofer, A. (2010). Comparing brain networks of different size and connectivity density using graph theory. *PLoS ONE* 5, e13701.** ◐
  - 질문: 크기와 밀도가 다른 뇌 네트워크를 공정하게 비교할 수 있는가?
  - 핵심: 여러 보정법을 검토했지만 어떤 방법도 완전히 편향 없는 비교를 보장하지 못했습니다.
  - 우리 논문에서: 널 정규화만으로는 부족해 노드 수·밀도를 직접 맞춘 이유입니다(Introduction, 3.2절).
  - 읽을 곳: 결론 절

- [ ] **Zalesky, A. et al. (2010). Whole-brain anatomical networks: does the choice of nodes matter? *NeuroImage* 50, 970–983.** ◐
  - 질문: 노드(뇌 영역)를 어떻게 나누느냐가 네트워크 지표를 바꾸는가?
  - 핵심: 같은 뇌라도 노드 수에 따라 small-world 지수 같은 값이 크게 바뀌었지만, "small-world인가" 같은 이분법적 결론은 유지됐습니다.
  - 우리 논문에서: 방향의 보존과 값의 보존을 나눠 봐야 한다는 근거입니다(Introduction).
  - 읽을 곳: 노드 수에 따른 지표 변화 그림

- [ ] **Rubinov, M. & Sporns, O. (2010). Complex network measures of brain connectivity: uses and interpretations. *NeuroImage* 52, 1059–1069.** ◐
  - 질문: 뇌 네트워크 지표를 어떻게 정의하고 해석해야 하는가?
  - 핵심: 클러스터링, 모듈성, 중심성, matching index 등의 표준 정의를 정리했습니다(Brain Connectivity Toolbox의 기반).
  - 우리 논문에서: 생성 모델의 matching index와 매개 중심성 정의의 출처입니다(3.6절).
  - 읽을 곳: 필요한 지표의 정의 표만

- [ ] **Murphy, A.H. (1988). Skill scores based on the mean square error and their relationships to the correlation coefficient. *Mon Weather Rev* 116, 2417–2424.** `【확인 필요】`
  - 질문: 예보가 기준 예보(기후값)보다 얼마나 나은지를 어떻게 하나의 점수로 나타내는가?
  - 핵심: skill score = (기준 오차 − 예보 오차) / (기준 오차 − 완벽 예보 오차)의 성질을 정리했습니다.
  - 우리 논문에서: 전이 지수 TI가 같은 구조라는 근거입니다. TI의 기준은 ER, 완벽은 대상 자신의 규칙입니다(3.7절, `docs/rule_transfer_methods.md`).
  - 읽을 곳: skill score 정의 부분만

---

**아직 확인할 것 (투고 전)**
- 저자 확인: Akarca et al. (2023), "Coming up short"(2025), 유충-성체 비교(Netw Neurosci 2025)
- 권·번호 확인: Skuhersky et al. (2022), Murphy (1988)
- ◐ 표시 논문의 권·쪽 번호는 DOI나 출판사 페이지로 한 번씩 확인
