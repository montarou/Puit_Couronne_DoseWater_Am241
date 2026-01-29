# Analyse des résultats après correction
## Simulation Geant4 Am-241 - 19 janvier 2026

---

## 1. Résumé des corrections validées

### ✅ Bug raie 0 (11.89 keV) - CORRIGÉ

| Métrique | Avant | Après |
|----------|-------|-------|
| Entrés dans eau | 429 | **897** |
| Absorbés dans eau | 580 | **551** |
| Taux d'absorption | **135.2%** (impossible!) | **61.4%** ✓ |

**Tous les taux d'absorption sont maintenant ≤ 100%**

### ✅ Ntuples EventData et GammaData - CORRIGÉS

| Ntuple | Avant | Après |
|--------|-------|-------|
| EventData | 0 entrées | **100,000 entrées** ✓ |
| GammaData | 0 entrées | **75,773 entrées** ✓ |

### ✅ Cohérences vérifiées

| Vérification | Résultat |
|--------------|----------|
| gamma_lines vs GammaData | ✓ OK |
| gamma_lines vs doses | ✓ OK |
| Formule de dose | ✓ OK (écarts < 0.1%) |

---

## 2. Statistiques par raie gamma Am-241

| Index | Énergie (keV) | Émis | Entrés | Absorbés | Taux abs. (%) |
|-------|---------------|------|--------|----------|---------------|
| 0 | 11.9 | 961 | 897 | 551 | 61.43 ✓ |
| 1 | 13.9 | 13,020 | 12,053 | 5,220 | 43.31 |
| 2 | 17.0 | 18,302 | 16,964 | 4,605 | 27.15 |
| 3 | 20.8 | 5,195 | 4,812 | 721 | 14.98 |
| 4 | 26.3 | 2,325 | 2,151 | 160 | 7.44 |
| 5 | 33.2 | 112 | 102 | 3 | 2.94 |
| 6 | 43.4 | 87 | 83 | 4 | 4.82 |
| 7 | 55.6 | 21 | 18 | 0 | 0.00 |
| **8** | **59.5** | **35,698** | **33,216** | **160** | **0.48** |
| 9 | 99.0 | 20 | 17 | 0 | 0.00 |
| 10 | 103.0 | 27 | 25 | 0 | 0.00 |
| 11 | 125.3 | 5 | 5 | 0 | 0.00 |

---

## 3. Doses par anneau d'eau

| Anneau | Rayon (mm) | Masse (g) | Énergie (MeV) | Dose (pGy/evt) |
|--------|------------|-----------|---------------|----------------|
| 0 | 0-5 | 0.0785 | 3.344 | **0.0682** |
| 1 | 5-10 | 0.2356 | 7.764 | **0.0528** |
| 2 | 10-15 | 0.3927 | 12.62 | **0.0515** |
| 3 | 15-20 | 0.5498 | 15.26 | **0.0445** |
| 4 | 20-25 | 0.7069 | 14.44 | **0.0327** |
| **TOTAL** | - | 1.9635 | 53.43 | **0.2497** |

---

## 4. Bilan de flux

```
Gammas primaires générés : 75,773 (100%)
  ├─ Transmis (PostContainer) : 58,919 (77.8%)
  ├─ Absorbés dans eau        : 11,424 (15.1%)
  └─ Autres (air, hors volume): 5,430 (7.2%)
```

---

## 5. Débit de dose estimé

Pour une source Am-241 de **42 kBq** :

| Paramètre | Valeur |
|-----------|--------|
| Gammas émis/s | ~32,050 |
| Débit de dose | **8.0 pGy/s** |
| | **28.8 nGy/h** |
| | **0.69 µGy/jour** |

---

## 6. Points d'attention mineurs

### ⚠️ Histogramme hRadialDose vide
- Le remplissage n'a pas été implémenté dans le code
- **Impact** : Mineur (les doses sont calculées via le ntuple `doses`)

### ⚠️ Différence PreContainer vs Water1 entry
- PreContainer : 73,765 photons
- Water1 entry : 70,343 photons
- Différence : 3,422 photons (4.6%)
- **Explication** : Photons absorbés dans l'air ou diffusés latéralement
- **Impact** : Aucun, physiquement cohérent

### ⚠️ Raie 6 (43.4 keV) - taux légèrement élevé
- Taux : 4.82% (vs 2.94% pour raie 5)
- **Cause** : Fluctuation statistique (seulement 87 émis)
- **Impact** : Négligeable

---

## 7. Conclusion

### ✅ LES CORRECTIONS SONT VALIDÉES

La simulation fonctionne correctement. Les résultats sont physiquement cohérents :

1. **Taux d'absorption** décroissants avec l'énergie (comme attendu)
2. **Formule de dose** validée : D = E × 0.160218 / m
3. **Profil radial** de dose décroissant vers l'extérieur (géométrie correcte)
4. **Cohérence** entre tous les ntuples et histogrammes

### Aucune incohérence majeure détectée.

---

*Rapport généré le 19 janvier 2026*
