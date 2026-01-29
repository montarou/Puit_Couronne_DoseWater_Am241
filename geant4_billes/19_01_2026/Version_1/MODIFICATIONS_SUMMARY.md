# RÉSUMÉ DES MODIFICATIONS - Simulation Geant4 Am-241

**Date** : 19 janvier 2026  
**Projet** : Simulation Monte Carlo de l'irradiation d'eau par une source Am-241 (géométrie puits couronne)

---

## 1. CORRECTION MAJEURE : Eu-152 → Am-241

### Problème identifié
Le code contenait des références à l'Europium-152 alors que le spectre implémenté était celui de l'Américium-241.

### Fichiers modifiés

#### EventAction.hh / EventAction.cc
- `kNbGammaLines` : 13 → **12** (Am-241 a 12 raies, pas 13)
- Tableau `fGammaLineEnergies` : énergies Eu-152 remplacées par Am-241
- Tableau `fGammaLineNames` : noms des raies mis à jour

**Anciennes énergies (Eu-152)** :
```
39.52, 40.12, 121.78, 244.70, 344.28, 411.12, 443.97, 
778.90, 867.38, 964.08, 1085.87, 1112.07, 1408.01 keV
```

**Nouvelles énergies (Am-241)** :
```
11.89, 13.9, 17.0, 20.8, 26.3446, 33.1963, 43.420, 
55.56, 59.5409, 98.97, 102.98, 125.30 keV
```

#### RunAction.hh / RunAction.cc
- Commentaires : "Eu-152" → "Am-241"
- `fMeanGammasPerDecay` : 2.03 → **0.7631** (rendement gamma Am-241)
- `fSourcePosZ` : 75.0 mm → **73.5 mm**
- En-têtes des tableaux de statistiques mis à jour

#### DetectorConstruction.hh / DetectorConstruction.cc
- Commentaires de géométrie mis à jour
- Messages d'affichage corrigés

#### SteppingAction.cc
- Message de démarrage : "Eu-152" → "Am-241"

---

## 2. CORRECTION MAJEURE : Épaisseur des anneaux d'eau

### Problème identifié
Dans `RunAction.cc`, le calcul des masses utilisait une épaisseur de **5 mm** alors que les anneaux n'ont que **1 mm** d'épaisseur (z = 102-103 mm).

### Fichier modifié : RunAction.cc (ligne 274)

**Avant** :
```cpp
const G4double waterThickness = 5.0 * mm;
```

**Après** :
```cpp
const G4double waterThickness = 1.0 * mm;
```

### Impact sur les masses et doses

| Anneau | Rayon (mm) | Masse ancienne (g) | Masse corrigée (g) |
|--------|------------|--------------------|--------------------|
| 0      | 0-5        | 0.3927             | **0.0785**         |
| 1      | 5-10       | 1.1781             | **0.2356**         |
| 2      | 10-15      | 1.9635             | **0.3927**         |
| 3      | 15-20      | 2.7489             | **0.5498**         |
| 4      | 20-25      | 3.5343             | **0.7069**         |
| TOTAL  |            | 9.8175             | **1.9635**         |

**Conséquence** : Les doses étaient **5× sous-estimées** car D = E/m et m était 5× trop grand.

---

## 3. CORRECTION : Comptage PreContainer

### Problème identifié
Le comptage des photons au plan PreContainer incluait tous les gammas (primaires + secondaires), créant une incohérence avec le comptage Water1 (primaires seulement).

### Fichier modifié : SteppingAction.cc (ligne 301)

**Avant** :
```cpp
if (particleName == "gamma" && pz > 0) {
```

**Après** :
```cpp
if (particleName == "gamma" && parentID == 0 && pz > 0) {
```

---

## 4. CORRECTION : Commentaires erronés

### SteppingAction.cc

**PreContainerPlane** (ligne 295) :
- Avant : "PLAN PRE-CONTAINER (dans eau, entrée)"
- Après : "PLAN PRE-CONTAINER (AIR, z = 99-100 mm, AVANT la surface de l'eau)"

**PostContainerPlane** (ligne 334) :
- Avant : "PLAN POST-CONTAINER (après eau, matériau eau)"
- Après : "PLAN POST-CONTAINER (POLYSTYRÈNE, z = 103-104 mm, APRÈS l'eau)"

---

## 5. TABLEAU DES RAIES Am-241 IMPLÉMENTÉES

| Index | Énergie (keV) | Nom            | Intensité (%) |
|-------|---------------|----------------|---------------|
| 0     | 11.89         | X_Ll (Np)      | 1.0           |
| 1     | 13.9          | X_Lα (Np)      | 13.0          |
| 2     | 17.0          | X_Lβ (Np)      | 18.5          |
| 3     | 20.8          | X_Lγ (Np)      | 5.16          |
| 4     | 26.3446       | γ 2,1          | 2.31          |
| 5     | 33.1963       | γ 1,0          | 0.12          |
| 6     | 43.420        | γ 4,2          | 0.07          |
| 7     | 55.56         | γ 6,4          | 0.02          |
| **8** | **59.5409**   | **γ 2,0**      | **35.92**     |
| 9     | 98.97         | γ 6,2          | 0.02          |
| 10    | 102.98        | γ 4,0          | 0.02          |
| 11    | 125.30        | γ 6,1          | 0.004         |

---

## 6. GÉOMÉTRIE DU DÉTECTEUR (rappel)

```
Source Am-241 (z = 73.5 mm)
        │
        ▼ gammas (+z)
┌───────────────────────────┐
│  PreContainerPlane        │  z = 99-100 mm   │ AIR
├───────────────────────────┤
│  Water1 (uniforme)        │  z = 100-102 mm  │ EAU (2 mm)
├───────────────────────────┤
│  WaterRing_0..4 (anneaux) │  z = 102-103 mm  │ EAU (1 mm) ← MESURE DOSE
├───────────────────────────┤
│  PostContainerPlane       │  z = 103-104 mm  │ POLYSTYRÈNE
├───────────────────────────┤
│  TungstenFoil             │  z = 104-104.05  │ TUNGSTÈNE (50 µm)
└───────────────────────────┘
```

---

## 7. FORMULE DE DOSE

```
D (nGy) = E (MeV) × 0.160218 / m (g)
```

où `0.160218 = 1.60218×10⁻¹³ J/MeV × 10³ g/kg × 10⁹ nGy/Gy`

---

## 8. FICHIERS INCLUS DANS L'ARCHIVE

### Sources C++ (16 fichiers)
- ActionInitialization.cc / .hh
- DetectorConstruction.cc / .hh
- EventAction.cc / .hh
- Logger.cc / .hh
- PhysicsList.cc / .hh
- PrimaryGeneratorAction.cc / .hh
- RunAction.cc / .hh
- SteppingAction.cc / .hh

### Documentation
- MODIFICATIONS_SUMMARY.md (ce fichier)
- ROOT_STRUCTURE.md (description histogrammes/ntuples)
- ROOT_STRUCTURE.tex / .pdf (version LaTeX)

---

## 9. RECOMPILATION

Après remplacement des fichiers sources :

```bash
cd build
cmake ..
make -j$(nproc)
./simulation run.mac
```

---

## 10. RÉSULTATS ATTENDUS APRÈS CORRECTION

1. **Tableau des raies gamma** : sera correctement rempli avec les statistiques Am-241
2. **Doses** : seront **5× plus élevées** qu'avant (valeurs correctes)
3. **Comptages PreContainer** : cohérents avec Water1 (primaires seulement)
4. **Messages de log** : référenceront Am-241 au lieu de Eu-152

---

*Document généré automatiquement - Simulation Geant4 Am-241*
