from enum import IntFlag
class Propagation(IntFlag):
    NONE=0,
    TRANS_SYSTEM_DEFAULT=1,
    TRANS_LANGEVIN=2,
    TRANS_VS_RELATIVE=4,
    TRANS_LB_MOMENTUM_EXCHANGE=8,
    TRANS_LB_TRACER=16,
    TRANS_BROWNIAN=32,
    TRANS_STOKESIAN=64,
    ROT_LANGEVIN=128,
    ROT_VS_RELATIVE=256,
    ROT_BROWNIAN=512
