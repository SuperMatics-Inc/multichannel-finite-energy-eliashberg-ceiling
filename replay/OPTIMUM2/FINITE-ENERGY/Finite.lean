import Std

namespace FiniteEnergyOptimum

theorem head_square (g c x y : Int) :
    g * (g * x ^ 2 - 2 * c * x * y + g * y ^ 2) =
      (g * x - c * y) ^ 2 + (g ^ 2 - c ^ 2) * y ^ 2 := by
  grind

theorem decoupled_head_criterion (px py r t : Int) :
    (2 * t * (1 + t) * px + r) * (2 * t * (1 + t) * r * py + 1) - r * (1 + 2 * t) ^ 2 =
      2 * t * (1 + t) * (px + r ^ 2 * py + 2 * t * (1 + t) * r * px * py - 2 * r) := by
  grind

theorem reflected_tail_decomposition (a b x y : Int) :
    (a - b) * (x ^ 2 + y ^ 2) - 2 * (a + b) * x * y =
      (a - b) * (x - y) ^ 2 - 4 * b * x * y := by
  grind

theorem cascade01_width_diff (m : Int) :
    ((m + 1) ^ 2 - m ^ 2 = 2 * m + 1) ∧
    ((m + 2) ^ 2 - (m - 1) ^ 2 = 3 * (2 * m + 1)) := by
  grind

theorem stage0_transfer_identity (t m : Int) :
    (2 * m + 1) * ((2 * t + (m + 1) ^ 2 + m ^ 2) ^ 2 - ((m + 1) ^ 2 - m ^ 2) ^ 2) =
      4 * (t + m ^ 2) * (t + (m + 1) ^ 2) * ((m + 1) ^ 2 - m ^ 2) := by
  grind

theorem stage1_transfer_identity (t k : Int) :
    3 * (2 * k + 1) * ((2 * t + (k + 2) ^ 2 + (k - 1) ^ 2) ^ 2 - ((k + 2) ^ 2 - (k - 1) ^ 2) ^ 2) =
      4 * (t + (k - 1) ^ 2) * (t + (k + 2) ^ 2) * ((k + 2) ^ 2 - (k - 1) ^ 2) := by
  grind
theorem stage0_completion_numerator (wra p : Int) :
    2 * wra * p + p ^ 2 - (p ^ 2 + 4 * p) =
      2 * p * (wra - 2) := by
  grind

theorem stage1_completion_numerator (wrb p : Int) :
    3 * (2 * wrb * p + p ^ 2) - (3 * p ^ 2 + 4 * p) =
      2 * p * (3 * wrb - 2) := by
  grind

theorem exact_budgets_and_ceilings :
    ((2 : Rat) / (3 ^ 2) = 2 / 9) ∧
    ((1 : Rat) / (3 ^ 2) = 1 / 9) ∧
    ((7 : Rat) / 10 * (1 / 6) + (3 / 10) * (1 / 3) < 7 / 32) ∧
    ((1 : Rat) / 12 + 1 / 100 = 7 / 75) ∧
    ((9 : Rat) / 17 - (7 / 75) / (617 / 1000) = 11899 / 31467) ∧
    ((11899 : Rat) / 31467 > 317 / 1000) ∧
    ((2 : Rat) / 5 + (3 / 4) * (617 / 8000 - 1 / 9) / (617 / 1000) = 106229 / 296160) ∧
    ((106229 : Rat) / 296160 > 317 / 1000) ∧
    ((1 : Rat) / 3 + (4 / 9) * (617 / 6000 - 1 / 9) / (617 / 1000) = 16361 / 49977) ∧
    ((16361 : Rat) / 49977 > 317 / 1000) ∧
    ((19 : Rat) / 64 < 119 / 400) ∧
    ((19 : Rat) / 60 < 317 / 1000) ∧
    (2 * ((1 : Rat) - 19 / 64) = 45 / 32) ∧
    (2 * ((1 : Rat) - 19 / 60) = 41 / 30) ∧
    ((41 : Rat) / 30 < 45 / 32) ∧
    ((45 : Rat) / 32 < 377 / 250) ∧
    ((377 : Rat) / 250 < 413 / 250) ∧
    ((1 : Rat) - (1 - 119 / 400) / (1 - 19 / 64) = 1 / 1125) ∧
    ((1 : Rat) - (1 - 317 / 1000) / (1 - 19 / 60) = 1 / 2050) := by
  native_decide

end FiniteEnergyOptimum
