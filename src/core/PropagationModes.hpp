#ifndef PROPAGATIONMODES_HPP
#define PROPAGATIONMODES_HPP

enum PropagationMode {
  NONE = 0,
  TRANS_SYSTEM_DEFAULT = 1,
  TRANS_LANGEVIN = 2,
  TRANS_VS_RELATIVE = 4,
  TRANS_LB_MOMENTUM_EXCHANGE = 8,
  TRANS_LB_TRACER = 16,
  TRANS_BROWNIAN = 32,
  TRANS_STOKESIAN = 64,
  ROT_LANGEVIN = 128,
  ROT_VS_RELATIVE = 256,
  ROT_BROWNIAN = 512
};

static bool is_valid_propagation_combination(int propagation) {
  if (propagation == 0)
    return true; // make sure if true or false
  if (!(propagation & (propagation - 1)))
    return true; // only one trans or rot
  if (propagation == TRANS_LANGEVIN + ROT_LANGEVIN ||
      propagation == TRANS_VS_RELATIVE + ROT_VS_RELATIVE ||
      propagation == TRANS_BROWNIAN + ROT_BROWNIAN)
    return true; // same mode for translation and rotation
  if (propagation == TRANS_VS_RELATIVE + ROT_LANGEVIN ||
      propagation == TRANS_LANGEVIN + ROT_VS_RELATIVE)
    return true;
  if (propagation == TRANS_LB_MOMENTUM_EXCHANGE + TRANS_VS_RELATIVE)
    return true;
  if (propagation ==
          TRANS_LB_MOMENTUM_EXCHANGE + TRANS_VS_RELATIVE + ROT_LANGEVIN ||
      propagation ==
          TRANS_LB_MOMENTUM_EXCHANGE + TRANS_VS_RELATIVE + ROT_VS_RELATIVE)
    return true;
  return false;
}

#endif
