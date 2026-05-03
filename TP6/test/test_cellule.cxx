/**
 * @file test_cellule.cxx
 * @brief Tests unitaires pour la classe cellule.
 *
 * Ce fichier contient un ensemble de tests unitaires basés sur GoogleTest
 * permettant de vérifier le bon fonctionnement de la classe cellule.
 *
 * Les tests couvrent :
 * - le constructeur par défaut ;
 * - le constructeur avec particules et voisins ;
 * - l'ajout de particules ;
 * - le vidage d'une cellule ;
 * - l'ajout de voisins ;
 * - l'accès aux particules et aux voisins.
 */

#include <gtest/gtest.h>
#include <vector>
#include "cellule.hxx"
#include "particule.hxx"
#include "vecteur.hxx"

/**
 * @brief Vérifie que le constructeur par défaut crée une cellule vide sans voisins.
 */
TEST(CelluleTest, ConstructeurParDefaut) {
    cellule c;

    EXPECT_TRUE(c.getParticules().empty());
    EXPECT_TRUE(c.getVoisins().empty());
    EXPECT_EQ(c.getParticules().size(), 0);
    EXPECT_EQ(c.getVoisins().size(), 0);
}

/**
 * @brief Vérifie que le constructeur avec paramètres initialise correctement
 * les particules et les voisins fournis.
 */
TEST(CelluleTest, ConstructeurAvecParametres) {
    particule p1(1, 0, 1.0, vecteur(1.0, 0.0, 0.0), vecteur());
    particule p2(2, 0, 2.0, vecteur(0.0, 1.0, 0.0), vecteur());

    cellule voisin1;
    cellule voisin2;

    std::vector<particule*> particules = {&p1, &p2};
    std::vector<cellule*> voisins = {&voisin1, &voisin2};

    cellule c(particules, voisins);

    ASSERT_EQ(c.getParticules().size(), 2);
    ASSERT_EQ(c.getVoisins().size(), 2);

    EXPECT_EQ(c.getParticules()[0], &p1);
    EXPECT_EQ(c.getParticules()[1], &p2);

    EXPECT_EQ(c.getVoisins()[0], &voisin1);
    EXPECT_EQ(c.getVoisins()[1], &voisin2);
}

/**
 * @brief Vérifie qu'une particule peut être ajoutée correctement dans une cellule.
 */
TEST(CelluleTest, AjouteParticule) {
    cellule c;
    particule p(1, 0, 1.0, vecteur(1.0, 2.0, 3.0), vecteur());

    c.ajoute_particule(&p);

    ASSERT_EQ(c.getParticules().size(), 1);
    EXPECT_EQ(c.getParticules()[0], &p);
}

/**
 * @brief Vérifie que plusieurs particules peuvent être ajoutées dans une cellule.
 */
TEST(CelluleTest, AjoutePlusieursParticules) {
    cellule c;
    particule p1(1, 0, 1.0, vecteur(), vecteur());
    particule p2(2, 0, 1.0, vecteur(), vecteur());

    c.ajoute_particule(&p1);
    c.ajoute_particule(&p2);

    ASSERT_EQ(c.getParticules().size(), 2);
    EXPECT_EQ(c.getParticules()[0], &p1);
    EXPECT_EQ(c.getParticules()[1], &p2);
}

/**
 * @brief Vérifie que la méthode vide supprime toutes les particules de la cellule.
 */
TEST(CelluleTest, Vide) {
    cellule c;
    particule p1(1, 0, 1.0, vecteur(), vecteur());
    particule p2(2, 0, 1.0, vecteur(), vecteur());

    c.ajoute_particule(&p1);
    c.ajoute_particule(&p2);

    ASSERT_EQ(c.getParticules().size(), 2);

    c.vide();

    EXPECT_TRUE(c.getParticules().empty());
    EXPECT_EQ(c.getParticules().size(), 0);
}

/**
 * @brief Vérifie qu'une cellule voisine peut être ajoutée correctement.
 */
TEST(CelluleTest, AjouteVoisin) {
    cellule c;
    cellule voisin;

    c.ajoute_voisin(&voisin);

    ASSERT_EQ(c.getVoisins().size(), 1);
    EXPECT_EQ(c.getVoisins()[0], &voisin);
}

/**
 * @brief Vérifie que plusieurs voisins peuvent être ajoutés correctement.
 */
TEST(CelluleTest, AjoutePlusieursVoisins) {
    cellule c;
    cellule voisin1;
    cellule voisin2;

    c.ajoute_voisin(&voisin1);
    c.ajoute_voisin(&voisin2);

    ASSERT_EQ(c.getVoisins().size(), 2);
    EXPECT_EQ(c.getVoisins()[0], &voisin1);
    EXPECT_EQ(c.getVoisins()[1], &voisin2);
}

/**
 * @brief Vérifie que le getter non constant permet bien d'accéder
 * au conteneur de particules.
 */
TEST(CelluleTest, GetterParticulesNonConst) {
    cellule c;
    particule p(1, 0, 1.0, vecteur(), vecteur());

    c.ajoute_particule(&p);

    ASSERT_EQ(c.getParticules().size(), 1);
    EXPECT_EQ(c.getParticules()[0], &p);
}

/**
 * @brief Vérifie que le getter non constant permet bien d'accéder
 * au conteneur de voisins.
 */
TEST(CelluleTest, GetterVoisinsNonConst) {
    cellule c;
    cellule voisin;

    c.getVoisins().push_back(&voisin);

    ASSERT_EQ(c.getVoisins().size(), 1);
    EXPECT_EQ(c.getVoisins()[0], &voisin);
}
