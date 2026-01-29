# CORRECTIONS APPLIQUÉES - 19 janvier 2026

## Problèmes identifiés et corrections

### 1. BUG CRITIQUE : Raie 0 - Taux d'absorption > 100%

**Symptôme** : Raie 0 (11.89 keV) montrait 580 absorbés pour seulement 429 entrés (taux = 135%)

**Cause** : Les gammas absorbés dans `Water1` (z=100-102mm) n'étaient pas toujours comptés comme "entrés" car `RecordContainerEntry()` pouvait ne pas être appelé dans certains cas de timing.

**Correction** (`EventAction.cc` - fonction `RecordGammaAbsorbed`) :
```cpp
if (volumeName.find("Water") != std::string::npos) {
    fPrimaryGammas[it->second].absorbedInWater = true;
    
    // CORRECTION: Si absorbé dans l'eau, alors forcément entré!
    fPrimaryGammas[it->second].enteredWater = true;  // ← AJOUTÉ
}
```

**Logique** : On ne peut pas être absorbé dans l'eau sans y être entré !

---

### 2. Ntuples EventData et GammaData vides

**Symptôme** : Les ntuples `EventData` (Ntuple 0) et `GammaData` (Ntuple 2) avaient 0 entrées.

**Cause** : Les fonctions de remplissage n'existaient pas.

**Corrections** :

#### `RunAction.hh` - Ajout des déclarations :
```cpp
void FillEventDataNtuple(G4int eventID, G4double edepTotal_keV,
                         const std::array<G4double, 5>& ringEdep_keV,
                         G4int nGammaEmitted, G4int nGammaWater);

void FillGammaDataNtuple(G4int eventID, G4double energy_keV,
                         G4int lineID, G4bool reachedWater, G4bool absorbed);
```

#### `RunAction.cc` - Ajout des implémentations :
```cpp
void RunAction::FillEventDataNtuple(...) {
    // Remplit le Ntuple 0 avec les données de l'événement
    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleDColumn(0, 1, edepTotal_keV);
    // ... colonnes 2-8
    analysisManager->AddNtupleRow(0);
}

void RunAction::FillGammaDataNtuple(...) {
    // Remplit le Ntuple 2 avec les données du gamma
    analysisManager->FillNtupleIColumn(2, 0, eventID);
    // ... colonnes 1-4
    analysisManager->AddNtupleRow(2);
}
```

#### `EventAction.cc` - Ajout des appels dans `EndOfEventAction()` :
```cpp
// Remplissage EventData (1 entrée par événement)
fRunAction->FillEventDataNtuple(eventID, totalDeposit/keV, ringEdep_keV,
                                 fPrimaryGammas.size(), fGammasEnteredWater.size());

// Remplissage GammaData (1 entrée par gamma primaire)
for (const auto& gamma : fPrimaryGammas) {
    fRunAction->FillGammaDataNtuple(eventID, gamma.energyInitial/keV,
                                     gamma.gammaLineIndex,
                                     gamma.enteredWater, gamma.absorbedInWater);
}
```

---

## Fichiers modifiés

| Fichier | Modifications |
|---------|---------------|
| `EventAction.cc` | `RecordGammaAbsorbed()` + appels ntuples |
| `RunAction.cc` | `FillEventDataNtuple()` + `FillGammaDataNtuple()` |
| `RunAction.hh` | Déclarations des 2 nouvelles fonctions |

---

## Résultats attendus après recompilation

1. **Raie 0** : Taux d'absorption ≤ 100% (corrigé)
2. **EventData** : 100,000 entrées (1 par événement)
3. **GammaData** : ~76,000 entrées (1 par gamma primaire)

---

## Recompilation

```bash
cd build
cmake ..
make -j$(nproc)
./simulation run.mac
```

---

*Corrections générées le 19 janvier 2026*
