from tf_dpmd_kit import plm

float_deltag, float_pka = plm.prob_to_pka(
    float_prob = 138540/3861460,
    float_T = 330,
    float_volume = 15.6793091675**3,
)
print(float_deltag, float_pka)
