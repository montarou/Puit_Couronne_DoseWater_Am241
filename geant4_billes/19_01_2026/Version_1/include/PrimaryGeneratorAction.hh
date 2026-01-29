#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

/// @brief Génération des particules primaires selon le spectre Am-241
///
/// Cette classe génère des photons (X et gamma) selon le spectre de l'Américium-241.
/// 
/// IMPORTANT - Philosophie de la simulation :
/// ------------------------------------------
/// L'Am-241 se désintègre à 100% par émission alpha vers Np-237.
/// Chaque désintégration émet :
///   - 1 particule alpha (NON SIMULÉE - mais prise en compte dans l'activité)
///   - Des électrons de conversion et Auger (NON SIMULÉS)
///   - Des photons X et gamma (SEULS SIMULÉS)
///
/// Le nombre moyen de photons X/gamma par désintégration est ~0.76
/// (somme des intensités des raies sélectionnées).
///
/// Un événement Geant4 = une désintégration Am-241
/// Mais certains événements peuvent ne générer aucun photon (événements "vides")
/// car les photons ne sont pas émis à chaque désintégration.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event*);
    
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // ═══════════════════════════════════════════════════════════════
    // ACCESSEURS POUR DIAGNOSTIC
    // ═══════════════════════════════════════════════════════════════
    G4int GetLastEventGammaCount() const { return fLastEventGammaCount; }

    // Accès au spectre (pour vérification)
    const std::vector<G4double>& GetGammaEnergies() const { return fGammaEnergies; }
    const std::vector<G4double>& GetGammaProbabilities() const { return fGammaProbabilities; }
    
    // ═══════════════════════════════════════════════════════════════
    // ACCESSEURS ET MODIFICATEURS POUR LE CÔNE D'ÉMISSION
    // ═══════════════════════════════════════════════════════════════
    void SetConeAngle(G4double angle) { fConeAngle = angle; }
    G4double GetConeAngle() const { return fConeAngle; }
    
    /// @brief Retourne le nombre moyen de photons X/gamma par désintégration Am-241
    /// Calculé à partir des raies sélectionnées (somme des intensités)
    /// Note: Chaque désintégration émet aussi 1 alpha (non simulé)
    static G4double GetMeanPhotonsPerDecay() { return 0.7631; }
    
    /// @brief Retourne le nombre de raies X/gamma implémentées
    static G4int GetNumberOfLines() { return 12; }
    
private:
    G4ParticleGun* fParticleGun;

    // ═══════════════════════════════════════════════════════════════
    // SPECTRE X/GAMMA Américium-241
    // Source: LNHB (Laboratoire National Henri Becquerel)
    // ═══════════════════════════════════════════════════════════════
    std::vector<G4double> fGammaEnergies;        // Énergies des raies (keV)
    std::vector<G4double> fGammaIntensities;     // Intensités relatives (%)
    std::vector<G4double> fGammaProbabilities;   // Probabilités d'émission (fraction 0-1)
    std::vector<G4String> fGammaNames;           // Noms des raies pour identification

    // ═══════════════════════════════════════════════════════════════
    // COMPTEUR POUR LE DERNIER ÉVÉNEMENT
    // ═══════════════════════════════════════════════════════════════
    G4int fLastEventGammaCount;

    // ═══════════════════════════════════════════════════════════════
    // PARAMÈTRES DE LA SOURCE
    // ═══════════════════════════════════════════════════════════════
    G4double fConeAngle;            // Demi-angle du cône d'émission
    G4ThreeVector fSourcePosition;  // Position de la source

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODE POUR GÉNÉRER UNE DIRECTION DANS LE CÔNE
    // ═══════════════════════════════════════════════════════════════
    void GenerateDirectionInCone(G4double coneAngle,
                                 G4double& theta,
                                 G4double& phi,
                                 G4ThreeVector& direction);
};

#endif
