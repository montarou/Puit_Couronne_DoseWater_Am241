#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fLastEventGammaCount(0),
  fConeAngle(45.*deg),
  fSourcePosition(0., 0., 73.5*mm)  // Position source (à ajuster selon géométrie)
{
    // Créer le particle gun
    fParticleGun = new G4ParticleGun(1);
    
    // Configuration par défaut : gamma
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(fSourcePosition);
    
    // ═══════════════════════════════════════════════════════════════
    // SPECTRE X/GAMMA AMÉRICIUM-241
    // ═══════════════════════════════════════════════════════════════
    // Source: LNHB (Laboratoire National Henri Becquerel)
    // http://www.lnhb.fr/nuclides/Am-241_tables.pdf
    //
    // L'Am-241 décroît à 100% par émission alpha vers Np-237.
    // Les photons X proviennent de la réorganisation du cortège 
    // électronique du Np-237 (raies X L du neptunium).
    // Les photons gamma proviennent de la désexcitation des niveaux
    // excités du Np-237.
    //
    // NOTE IMPORTANTE:
    // - Les alpha (100%) ne sont PAS simulés
    // - Les électrons de conversion et Auger ne sont PAS simulés
    // - Seuls les photons X et gamma sont simulés
    // - L'activité de 42 kBq correspond à 42000 désintégrations/s
    //   (chacune émettant 1 alpha, même si on ne le simule pas)
    // ═══════════════════════════════════════════════════════════════
    
    // ───────────────────────────────────────────────────────────────
    // RAIES X L DU NEPTUNIUM (réorganisation électronique après α)
    // Intensité totale X_L = 37.66% répartie en sous-groupes
    // ───────────────────────────────────────────────────────────────
    // Les raies X L sont regroupées par sous-groupes avec énergies moyennes
    
    fGammaEnergies.push_back(11.89);   // X_Ll (Np)
    fGammaIntensities.push_back(1.0);
    fGammaNames.push_back("X_Ll");
    
    fGammaEnergies.push_back(13.9);    // X_Lalpha (Np) - moyenne Lalpha1,2
    fGammaIntensities.push_back(13.0);
    fGammaNames.push_back("X_Lalpha");
    
    fGammaEnergies.push_back(17.0);    // X_Lbeta (Np) - moyenne des Lbeta
    fGammaIntensities.push_back(18.5);
    fGammaNames.push_back("X_Lbeta");
    
    fGammaEnergies.push_back(20.8);    // X_Lgamma (Np) - moyenne des Lgamma
    fGammaIntensities.push_back(5.16);
    fGammaNames.push_back("X_Lgamma");
    
    // ───────────────────────────────────────────────────────────────
    // RAIES GAMMA (désexcitation du Np-237)
    // Sélection des raies avec intensité > 0.01%
    // ───────────────────────────────────────────────────────────────
    
    fGammaEnergies.push_back(26.3446); // gamma 2,1 (Np)
    fGammaIntensities.push_back(2.31);
    fGammaNames.push_back("gamma_26keV");
    
    fGammaEnergies.push_back(33.1963); // gamma 1,0 (Np)
    fGammaIntensities.push_back(0.1215);
    fGammaNames.push_back("gamma_33keV");
    
    fGammaEnergies.push_back(43.420);  // gamma 4,2 (Np)
    fGammaIntensities.push_back(0.0669);
    fGammaNames.push_back("gamma_43keV");
    
    fGammaEnergies.push_back(55.56);   // gamma 6,4 (Np)
    fGammaIntensities.push_back(0.0181);
    fGammaNames.push_back("gamma_56keV");
    
    fGammaEnergies.push_back(59.5409); // gamma 2,0 (Np) - RAIE PRINCIPALE
    fGammaIntensities.push_back(35.92);
    fGammaNames.push_back("gamma_59keV");
    
    fGammaEnergies.push_back(98.97);   // gamma 6,2 (Np)
    fGammaIntensities.push_back(0.0203);
    fGammaNames.push_back("gamma_99keV");
    
    fGammaEnergies.push_back(102.98);  // gamma 4,0 (Np)
    fGammaIntensities.push_back(0.0195);
    fGammaNames.push_back("gamma_103keV");
    
    fGammaEnergies.push_back(125.30);  // gamma 6,1 (Np)
    fGammaIntensities.push_back(0.0041);
    fGammaNames.push_back("gamma_125keV");
    
    // ═══════════════════════════════════════════════════════════════
    // CALCUL DES PROBABILITÉS
    // ═══════════════════════════════════════════════════════════════
    // Pour chaque raie, probabilité = intensité / 100
    // (car on génère tous les photons indépendamment)
    
    G4double totalIntensity = 0.;
    for (const auto& intensity : fGammaIntensities) {
        totalIntensity += intensity;
    }
    
    fGammaProbabilities.resize(fGammaEnergies.size());
    for (size_t i = 0; i < fGammaIntensities.size(); ++i) {
        fGammaProbabilities[i] = fGammaIntensities[i] / 100.;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // AFFICHAGE DU SPECTRE
    // ═══════════════════════════════════════════════════════════════
    
    G4cout << "\n╔═══════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║      PrimaryGeneratorAction: Spectre Am-241 initialisé           ║" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Source: LNHB (Laboratoire National Henri Becquerel)             ║" << G4endl;
    G4cout << "║  Désintégration: 100% alpha vers Np-237                          ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "║  *** ATTENTION: Seuls les photons X/gamma sont simulés ***       ║" << G4endl;
    G4cout << "║  Les alpha et électrons ne sont PAS générés mais sont            ║" << G4endl;
    G4cout << "║  pris en compte dans l'activité (42 kBq = 42000 désint./s)       ║" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Raies implémentées: " << fGammaEnergies.size() << "                                          ║" << G4endl;
    G4cout << "║  Intensité totale X+gamma: " << totalIntensity << "%                            ║" << G4endl;
    G4cout << "║  Photons moyens/désintégration: " << totalIntensity/100. << "                       ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "║  Raie principale: 59.54 keV (35.92%)                             ║" << G4endl;
    G4cout << "║  Raies X L (Np): 11.9-20.8 keV (37.66% total)                    ║" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Angle du cône: " << fConeAngle/deg << " degrés                                      ║" << G4endl;
    G4cout << "║  Position source: z = " << fSourcePosition.z()/mm << " mm                            ║" << G4endl;
    G4cout << "╚═══════════════════════════════════════════════════════════════════╝" << G4endl;
    
    // Détail des raies
    G4cout << "\n  Détail des raies:" << G4endl;
    G4cout << "  ─────────────────────────────────────────────────" << G4endl;
    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {
        G4cout << "    " << fGammaNames[i] << ": " 
               << fGammaEnergies[i] << " keV, I = " 
               << fGammaIntensities[i] << "%" << G4endl;
    }
    G4cout << "  ─────────────────────────────────────────────────\n" << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // ═══════════════════════════════════════════════════════════════
    // GÉNÉRATION DES PHOTONS X/GAMMA POUR UNE DÉSINTÉGRATION Am-241
    // ═══════════════════════════════════════════════════════════════
    // 
    // Chaque événement = une désintégration Am-241
    // Pour chaque raie, on tire si elle est émise (Bernoulli indépendant)
    // 
    // Note: La plupart des événements (~24%) ne génèreront aucun photon
    // car l'intensité totale n'est que de ~76%
    // ═══════════════════════════════════════════════════════════════
    
    fLastEventGammaCount = 0;
    
    // Pour chaque raie X/gamma, tirer si elle est émise
    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {
        G4double random = G4UniformRand();
        
        if (random < fGammaProbabilities[i]) {
            // Cette raie est émise
            
            // Énergie de la raie
            G4double energy = fGammaEnergies[i] * keV;
            
            // Générer une direction dans le cône
            G4double theta, phi;
            G4ThreeVector direction;
            GenerateDirectionInCone(fConeAngle, theta, phi, direction);
            
            // Configurer et tirer
            fParticleGun->SetParticleEnergy(energy);
            fParticleGun->SetParticleMomentumDirection(direction);
            fParticleGun->SetParticlePosition(fSourcePosition);
            fParticleGun->GeneratePrimaryVertex(anEvent);
            
            fLastEventGammaCount++;
        }
    }
    
    // Note: Si aucun photon n'a été émis, l'événement est "vide" côté photons
    // mais représente quand même une désintégration (avec émission alpha)
    // C'est physiquement correct: ~24% des désintégrations Am-241 n'émettent
    // pas de photon X/gamma détectable dans notre liste de raies
}

void PrimaryGeneratorAction::GenerateDirectionInCone(G4double coneAngle,
                                                     G4double& theta,
                                                     G4double& phi,
                                                     G4ThreeVector& direction)
{
    // Distribution uniforme sur la surface de la calotte sphérique
    // cos(theta) est uniforme entre cos(coneAngle) et 1
    G4double cosTheta = 1. - G4UniformRand() * (1. - std::cos(coneAngle));
    theta = std::acos(cosTheta);
    
    // phi est uniforme entre 0 et 2π
    phi = G4UniformRand() * 2. * CLHEP::pi;
    
    // Calculer la direction
    G4double sinTheta = std::sin(theta);
    direction.set(sinTheta * std::cos(phi),
                  sinTheta * std::sin(phi),
                  cosTheta);
}
