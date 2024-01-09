#ifndef M_MATRICES_H_INCLUDED
#define M_MATRICES_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>

class M_matrices{

public:
    M_matrices();
    M_matrices(int nb_lignes, int nb_col);
    M_matrices(int nb_lignes, int nb_col, std::vector<double> matrice_vector);

    std::vector<int> donner_index_non_nul();

    int get_nb_lignes();
    int get_nb_col();
    std::vector<double> get_matrice_vect();
    std::vector<double> get_matrice_inverse();
    bool is_inversible() const;


    void set_matrice_vect(std::vector<double> matrice_vector);
    void swap_lignes();
    void reduire();
    void echelonner();
    void echelonner_reduire();
    void initialise_matrice_inverse();





    void affiche(std::ostream& flux) const;
    bool estEgal(M_matrices const& autre) const;
    void recupere();
    M_matrices get_inverse();
    void set_nb_lignes_col(int nb_lignes, int nb_col);



    M_matrices& operator+=(M_matrices const& autre);
    M_matrices& operator-=(M_matrices const& autre);
    M_matrices& operator*=(double const& lambda);
    M_matrices& operator/=(double const& lambda);
    M_matrices operator*(M_matrices const& autre);




private:

    int m_nb_col;
    int m_nb_lignes;
    std::vector <double> m_matrice;
    std::vector <double> m_matrice_inverse;

};
    std::ostream& operator<<(std::ostream& flux, M_matrices const& matrice);
    std::istream& operator>>(std::istream& flux, M_matrices &matrice);
    void afficher_titre();


    M_matrices get_matrice_identite(int nb_lignes, int nb_col);
    //M_matrices get_inverse(M_matrices m);


    M_matrices operator+(M_matrices const& a, M_matrices const& b);
    M_matrices operator-(M_matrices const& a, M_matrices const& b);
    M_matrices operator*(M_matrices const& a, M_matrices const& b);
    M_matrices operator*(M_matrices const& a, double const& lambda);
    M_matrices operator/(M_matrices const& a, double const& lambda);

    bool operator==(M_matrices const& a, M_matrices const& b);
    bool operator!=(M_matrices const& a, M_matrices const& b);


#endif // M_MATRICES_H_INCLUDED
