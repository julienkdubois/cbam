CBAM -- Dynare models repository

Description des scripts et organisation

- `dynare/call_file/ss_residuals.m` : utilitaire utilisé pour ajuster les résidus du steady state afin de matcher le PIB.
- `dynare/cbam_estimation.mod` et `dynare/ss_pf.m` : fichiers dédiés à l'estimation bayésienne et au calcul du steady-state (estimation de paramètres).
- `dynare/pf_cbam.mod` : modèle en perfect-foresight pour simuler la transition CBAM ; fonctionne avec `ss_pf.m` pour résoudre l'état stationnaire.
- `dynare/compare_pf.m` : script de comparaison des trajectoires (sticky vs flexible) et de simulation pour différentes dates d'implémentation. Options ajoutées : horizon par défaut = 40, sauvegarde des figures, comparaison multi-dates.
- `dynare/old/` : versions historiques et expérimentales (exclues des commits récents par décision de maintenance).

Usage rapide

- Lancer une comparaison pour trois dates d'implémentation :
  ```matlab
  compare_pf([1,13,30])
  ```
- Pour utiliser les paramètres estimés et patcher `pf_cbam.mod`, répondre `y` à la question interactive dans `compare_pf`.

Notes

- Les modifications du 30/04/2026 ont harmonisé la gestion de `theta1` entre estimation et perfect-foresight (fix dans `ss_pf.m`) et adapté `compare_pf.m` pour des usages multi-dates et sauvegarde de graphiques.
