# Structure du fichier ROOT - Simulation Geant4 Am-241

**Fichier** : `output.root`  
**Date** : 19 janvier 2026

---

## 1. Histogrammes 1D (16 total)

| Index | Nom | Description | Range | Bins | Entrées |
|-------|-----|-------------|-------|------|---------|
| H1[0] | hGammaEmitted | Spectre γ émis | 0–2000 keV | 1000 | 75,773 |
| H1[1] | hGammaEnteringWater | Spectre γ entrant eau | 0–2000 keV | 1000 | 54,947 |
| H1[2] | hEdepWater | Dépôt énergie eau (step) | 0–250 keV | 500 | 12,783 |
| H1[3] | hEdepRing0 | Dépôt énergie anneau 0 | 0–200 keV | 200 | 751 |
| H1[4] | hEdepRing1 | Dépôt énergie anneau 1 | 0–200 keV | 200 | 1,929 |
| H1[5] | hEdepRing2 | Dépôt énergie anneau 2 | 0–200 keV | 200 | 3,034 |
| H1[6] | hEdepRing3 | Dépôt énergie anneau 3 | 0–200 keV | 200 | 3,590 |
| H1[7] | hEdepRing4 | Dépôt énergie anneau 4 | 0–200 keV | 200 | 3,479 |
| H1[8] | hRadialDose | Profil radial dose | 0–25 mm | 25 | **0** ⚠️ |
| H1[9] | hElectronSpectrum | Spectre e⁻ secondaires | 0–1500 keV | 500 | 8,655 |
| H1[10] | h_dose_ring0 | Dose/evt anneau 0 | variable | var. | 239 |
| H1[11] | h_dose_ring1 | Dose/evt anneau 1 | variable | var. | 609 |
| H1[12] | h_dose_ring2 | Dose/evt anneau 2 | variable | var. | 959 |
| H1[13] | h_dose_ring3 | Dose/evt anneau 3 | variable | var. | 1,135 |
| H1[14] | h_dose_ring4 | Dose/evt anneau 4 | variable | var. | 1,101 |
| H1[15] | h_dose_total | Dose totale/evt | 0–0.03 nGy | 500 | 3,973 |

### Conditions de remplissage

#### H1[0] – hGammaEmitted
- **Variable** : Énergie initiale du gamma (keV)
- **Condition** : Premier step d'un gamma primaire (`parentID==0 && stepNumber==1`)
- **Source** : `SteppingAction.cc` → `FillGammaEmittedSpectrum()`

#### H1[1] – hGammaEnteringWater
- **Variable** : Énergie du gamma entrant (keV)
- **Condition** : Gamma primaire entrant dans `WaterRing_*` pour la première fois
- **Source** : `SteppingAction.cc` → `FillGammaEnteringWater()`
- **Note** : Anti-double-comptage via `HasEnteredWater(trackID)`

#### H1[2] – hEdepWater
- **Variable** : Dépôt d'énergie par step (keV)
- **Condition** : `edep > 0` dans un volume `WaterRing_*`
- **Source** : `SteppingAction.cc` → `FillEdepWater()`

#### H1[3–7] – hEdepRing0..4
- **Variable** : Dépôt d'énergie par step (keV)
- **Condition** : `edep > 0` dans l'anneau i correspondant
- **Source** : `SteppingAction.cc` → `FillEdepRing()`

#### H1[8] – hRadialDose
- **⚠️ NON IMPLÉMENTÉ** – histogramme créé mais jamais rempli

#### H1[9] – hElectronSpectrum
- **Variable** : Énergie cinétique de l'électron (keV)
- **Condition** : Électron secondaire (`parentID != 0`) dans `WaterRing_*`
- **Source** : `SteppingAction.cc` → `FillElectronSpectrum()`

#### H1[10–14] – h_dose_ring0..4
- **Variable** : Dose par événement (nGy)
- **Condition** : Fin d'événement si `edep > 0` dans l'anneau
- **Formule** : `D = E × 0.160218 / m`
- **Source** : `RunAction.cc`

#### H1[15] – h_dose_total
- **Variable** : Dose totale par événement (nGy)
- **Condition** : Fin d'événement
- **Formule** : `D_tot = Σ D_i`

---

## 2. Histogrammes 2D (2 total)

| Index | Nom | Description | Range | Bins | Entrées |
|-------|-----|-------------|-------|------|---------|
| H2[0] | hEdepXY | Carte XY des dépôts | [-30,30]×[-30,30] mm | 100×100 | 53,427 |
| H2[1] | hEdepRZ | Carte RZ des dépôts | [0,30]×[98,108] mm | 50×50 | 53,427 |

### Conditions de remplissage

#### H2[0] – hEdepXY
- **Variables** : X (mm), Y (mm)
- **Poids** : `edep` (keV)
- **Condition** : `edep > 0` dans `WaterRing_*`

#### H2[1] – hEdepRZ
- **Variables** : r = √(x²+y²) (mm), Z (mm)
- **Poids** : `edep` (keV)
- **Condition** : `edep > 0` dans `WaterRing_*`

---

## 3. Ntuples (7 total)

### Ntuple 0 – EventData (100,000 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | EventID | Int | Numéro de l'événement |
| 1 | EdepTotal | Double | Énergie totale déposée (keV) |
| 2 | EdepRing0 | Double | Énergie déposée anneau 0 (keV) |
| 3 | EdepRing1 | Double | Énergie déposée anneau 1 (keV) |
| 4 | EdepRing2 | Double | Énergie déposée anneau 2 (keV) |
| 5 | EdepRing3 | Double | Énergie déposée anneau 3 (keV) |
| 6 | EdepRing4 | Double | Énergie déposée anneau 4 (keV) |
| 7 | NGammaEmitted | Int | Nombre de gammas primaires émis |
| 8 | NGammaWater | Int | Nombre de gammas entrés dans l'eau |

**Fréquence** : 1 entrée/événement  
**Remplissage** : `EndOfEventAction()` → `FillEventDataNtuple()`

---

### Ntuple 1 – StepData (12,783 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | EventID | Int | Numéro de l'événement |
| 1 | X | Double | Position X (mm) |
| 2 | Y | Double | Position Y (mm) |
| 3 | Z | Double | Position Z (mm) |
| 4 | Edep | Double | Dépôt d'énergie (keV) |
| 5 | RingID | Int | Index de l'anneau (0–4) |
| 6 | ParticleName | String | Nom de la particule |
| 7 | ProcessName | String | Nom du processus physique |

**Fréquence** : N entrées/événement (selon dépôts)  
**Condition** : `edep > 0` dans `WaterRing_*`  
**Remplissage** : `SteppingAction::UserSteppingAction()` → `FillStepNtuple()`

---

### Ntuple 2 – GammaData (75,773 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | EventID | Int | Numéro de l'événement |
| 1 | Energy | Double | Énergie initiale (keV) |
| 2 | LineID | Int | Index de la raie (0–11, -1 si inconnu) |
| 3 | ReachedWater | Int | 1 si entré dans l'eau, 0 sinon |
| 4 | Absorbed | Int | 1 si absorbé dans l'eau, 0 sinon |

**Fréquence** : 1 entrée/gamma primaire  
**Remplissage** : `EndOfEventAction()` → `FillGammaDataNtuple()`  
**Note** : `ReachedWater=1` inclut Water1 ET WaterRings (après correction)

---

### Ntuple 3 – gamma_lines (12 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | lineIndex | Int | Index de la raie (0–11) |
| 1 | energy_keV | Double | Énergie de la raie (keV) |
| 2 | emitted | Int | Nombre de gammas émis |
| 3 | enteredWater | Int | Nombre entrés dans l'eau |
| 4 | absorbedWater | Int | Nombre absorbés dans l'eau |
| 5 | waterAbsRate | Double | Taux d'absorption (%) |
| 6 | waterEntryRate | Double | Taux d'entrée (%) |

**Fréquence** : 12 entrées (1 par raie Am-241)  
**Remplissage** : `EndOfRunAction()` – statistiques agrégées

#### Raies Am-241 implémentées

| Index | Nom | Énergie (keV) | Intensité (%) |
|-------|-----|---------------|---------------|
| 0 | X_Ll (Np) | 11.89 | 1.0 |
| 1 | X_Lα (Np) | 13.9 | 13.0 |
| 2 | X_Lβ (Np) | 17.0 | 18.5 |
| 3 | X_Lγ (Np) | 20.8 | 5.16 |
| 4 | γ 2,1 | 26.3446 | 2.31 |
| 5 | γ 1,0 | 33.1963 | 0.12 |
| 6 | γ 4,2 | 43.420 | 0.07 |
| 7 | γ 6,4 | 55.56 | 0.02 |
| **8** | **γ 2,0 (PRINCIPALE)** | **59.5409** | **35.92** |
| 9 | γ 6,2 | 98.97 | 0.02 |
| 10 | γ 4,0 | 102.98 | 0.02 |
| 11 | γ 6,1 | 125.30 | 0.004 |

---

### Ntuple 4 – precontainer (100,000 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | eventID | Int | Numéro de l'événement |
| 1 | nPhotons | Int | Nombre de photons primaires (pz > 0) |
| 2 | sumEPhotons_keV | Double | Énergie totale photons (keV) |
| 3 | nElectrons | Int | Nombre d'électrons (pz > 0) |
| 4 | sumEElectrons_keV | Double | Énergie totale électrons (keV) |

**Fréquence** : 1 entrée/événement  
**Plan** : `PreContainerPlane` – AIR, z = 99–100 mm  
**Condition** : Particule entrant dans le plan avec pz > 0  
**Filtre photons** : `parentID == 0` (primaires seulement)

---

### Ntuple 5 – postcontainer (100,000 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | eventID | Int | Numéro de l'événement |
| 1 | nPhotons_fwd | Int | Photons transmis (pz > 0) |
| 2 | sumEPhotons_fwd_keV | Double | Énergie photons transmis (keV) |
| 3 | nPhotons_back | Int | Photons rétrodiffusés (pz < 0) |
| 4 | sumEPhotons_back_keV | Double | Énergie photons rétrodiff. (keV) |
| 5 | nElectrons_fwd | Int | Électrons transmis (pz > 0) |
| 6 | sumEElectrons_fwd_keV | Double | Énergie électrons transmis (keV) |
| 7 | nElectrons_back | Int | Électrons rétrodiffusés (pz < 0) |
| 8 | sumEElectrons_back_keV | Double | Énergie électrons rétrodiff. (keV) |

**Fréquence** : 1 entrée/événement  
**Plan** : `PostContainerPlane` – POLYSTYRÈNE, z = 103–104 mm  
**Filtre photons** : TOUS (primaires + Compton diffusés)

---

### Ntuple 6 – doses (100,000 entrées)

| Col | Nom | Type | Description |
|-----|-----|------|-------------|
| 0 | eventID | Int | Numéro de l'événement |
| 1 | dose_nGy_ring0 | Double | Dose anneau 0 (nGy) |
| 2 | dose_nGy_ring1 | Double | Dose anneau 1 (nGy) |
| 3 | dose_nGy_ring2 | Double | Dose anneau 2 (nGy) |
| 4 | dose_nGy_ring3 | Double | Dose anneau 3 (nGy) |
| 5 | dose_nGy_ring4 | Double | Dose anneau 4 (nGy) |
| 6 | dose_nGy_total | Double | Dose totale (nGy) |
| 7 | edep_keV_total | Double | Énergie totale déposée (keV) |
| 8 | nPrimaries | Int | Nombre de gammas primaires |
| 9 | nTransmitted | Int | Nombre de gammas transmis |
| 10 | nAbsorbed | Int | Nombre de gammas absorbés |

**Fréquence** : 1 entrée/événement  
**Remplissage** : `EndOfEventAction()` → `FillDosesNtuple()`  
**Formule dose** : `D (nGy) = E (MeV) × 0.160218 / m (g)`

---

## 4. Schéma de remplissage

| Moment | Fonction | Objets remplis |
|--------|----------|----------------|
| **Chaque step** | `UserSteppingAction()` | H1[0–9], H2[0–1], Ntuple 1, compteurs |
| **Fin événement** | `EndOfEventAction()` | H1[10–15], Ntuples 0, 2, 4, 5, 6 |
| **Fin run** | `EndOfRunAction()` | Ntuple 3 (gamma_lines) |

---

## 5. Géométrie de référence

```
Source Am-241 (z = 73.5 mm)
        │
        ▼ gammas (+z)
┌───────────────────────────────┐
│  PreContainerPlane            │  z = 99-100 mm   │ AIR
├───────────────────────────────┤
│  Water1 (uniforme)            │  z = 100-102 mm  │ EAU (2 mm)
├───────────────────────────────┤
│  WaterRing_0..4 (anneaux)     │  z = 102-103 mm  │ EAU (1 mm) ← MESURE DOSE
├───────────────────────────────┤
│  PostContainerPlane           │  z = 103-104 mm  │ POLYSTYRÈNE
├───────────────────────────────┤
│  TungstenFoil                 │  z = 104-104.05  │ W (50 µm)
└───────────────────────────────┘
```

### Dimensions des anneaux (épaisseur = 1 mm)

| Anneau | r_int (mm) | r_ext (mm) | Volume (cm³) | Masse (g) |
|--------|------------|------------|--------------|-----------|
| 0 | 0 | 5 | 0.0785 | 0.0785 |
| 1 | 5 | 10 | 0.2356 | 0.2356 |
| 2 | 10 | 15 | 0.3927 | 0.3927 |
| 3 | 15 | 20 | 0.5498 | 0.5498 |
| 4 | 20 | 25 | 0.7069 | 0.7069 |
| **TOTAL** | 0 | 25 | 1.9635 | 1.9635 |

---

*Document généré le 19 janvier 2026*
