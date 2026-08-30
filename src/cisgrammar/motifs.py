from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BASES = np.array(list("ACGT"))


@dataclass(frozen=True)
class Motif:
    name: str
    accession: str
    pwm: np.ndarray
    source: str
    source_url: str

    def __post_init__(self) -> None:
        pwm = np.asarray(self.pwm, dtype=np.float64)
        if pwm.ndim != 2 or pwm.shape[1] != 4:
            raise ValueError(f"{self.accession}: expected an L×4 matrix, received {pwm.shape}")
        if np.any(pwm < 0) or np.any(pwm.sum(axis=1) <= 0):
            raise ValueError(f"{self.accession}: matrix rows must be non-negative and non-empty")
        pwm = pwm / pwm.sum(axis=1, keepdims=True)
        object.__setattr__(self, "pwm", pwm)

    @property
    def length(self) -> int:
        return int(self.pwm.shape[0])

    def probabilities(self, temperature: float = 1.0) -> np.ndarray:
        if temperature <= 0:
            raise ValueError("motif temperature must be positive")
        scaled = np.power(np.clip(self.pwm, 1e-8, 1.0), 1.0 / temperature)
        return scaled / scaled.sum(axis=1, keepdims=True)

    def sample(self, rng: np.random.Generator, temperature: float = 1.0) -> str:
        matrix = self.probabilities(temperature)
        return "".join(rng.choice(BASES, p=row) for row in matrix)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def reverse_complement_pwm(pwm: np.ndarray) -> np.ndarray:
    matrix = np.asarray(pwm)
    return matrix[::-1, :][:, [3, 2, 1, 0]]


POU5F1 = Motif(
    name="POU5F1",
    accession="M05705_3.00",
    source="CIS-BP 3.00 (matrix supplied with the extended project)",
    source_url="https://cisbp.ccbr.utoronto.ca/",
    pwm=np.array(
        [
            [0.303787103377687, 0.231832139201638, 0.288024564994882, 0.176356192425793],
            [0.444977991196479, 0.129851940776311, 0.218554088301987, 0.206615979725223],
            [0.296343115124154, 0.0337697516930023, 0.0840632054176072, 0.585823927765237],
            [0.0207937618714386, 0.00249925022493252, 0.0, 0.976706987903629],
            [1.0, 0.0, 0.0, 0.0],
            [0.00010230179028133, 0.00010230179028133, 0.00030690537084399, 0.999488491048593],
            [0.000306873977086743, 0.000306873977086743, 0.999386252045827, 0.0],
            [0.00124797204542618, 0.60963434419069, 0.0207163359540746, 0.368401347809809],
            [0.464661406969099, 5.47885163269779e-05, 0.0, 0.535283804514574],
            [0.999590750971966, 0.000204624514016779, 0.00010231225700839, 0.00010231225700839],
            [0.999897656329956, 0.0, 0.000102343670044008, 0.0],
            [0.00428047289033836, 0.0, 0.0, 0.995719527109662],
            [0.00778785488958991, 0.0290812302839117, 0.441443217665615, 0.521687697160883],
            [0.649860316615671, 0.0527471065584675, 0.147399228415591, 0.14999334841027],
            [0.378266850068776, 0.199220541036222, 0.368332569157879, 0.0541800397371236],
            [0.0775844421699079, 0.239201637666325, 0.449846468781986, 0.233367451381781],
        ]
    ),
)


NANOG = Motif(
    name="NANOG",
    accession="M05219_3.00",
    source="CIS-BP 3.00 (matrix supplied with the extended project)",
    source_url="https://cisbp.ccbr.utoronto.ca/",
    pwm=np.array(
        [
            [0.167707791081034, 0.207444039562728, 0.291688356758633, 0.333159812597605],
            [0.0288970942366351, 0.23422700272917, 0.120404559319313, 0.616471343714882],
            [0.501370936153545, 0.0, 0.498629063846455, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.213439351418771, 0.187870283754287, 0.598690364826941],
            [0.0310373443983402, 0.331618257261411, 0.340082987551867, 0.297261410788382],
            [0.194831959416614, 0.153297400126823, 0.608750792644261, 0.0431198478123018],
            [0.113310758009898, 0.353737952591821, 0.285230528783537, 0.247720760614743],
        ]
    ),
)


MOTIFS = {
    "POU5F1_M05705_3.00": POU5F1,
    "NANOG_M05219_3.00": NANOG,
}


def get_motif(key: str) -> Motif:
    try:
        return MOTIFS[key]
    except KeyError as exc:
        raise KeyError(f"Unknown motif {key!r}; available motifs: {sorted(MOTIFS)}") from exc
