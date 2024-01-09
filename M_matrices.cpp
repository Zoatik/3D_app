#include <iostream>
#include "M_matrices.h"
#include <cmath>
#include <vector>
#include <stdlib.h>

using namespace std;


//constructeurs
    M_matrices::M_matrices(): m_nb_lignes(0),m_nb_col(0),m_matrice(0),m_matrice_inverse(0) ///on peut optimiser pour initialiser directement m_matrice_inverse
    {

    }
    M_matrices::M_matrices(int nb_lignes, int nb_col): m_nb_lignes(nb_lignes),m_nb_col(nb_col),m_matrice(0), m_matrice_inverse(0)
    {

    }
    M_matrices::M_matrices(int nb_lignes, int nb_col, vector<double> matrice_vector) : m_nb_lignes(nb_lignes), m_nb_col(nb_col),m_matrice(matrice_vector), m_matrice_inverse(0)
    {

    }


///fonctions hors classe


M_matrices get_matrice_identite(int nb_lignes, int nb_col)
{
    vector <double> matrice_vect;
    M_matrices matrice_identite(nb_lignes, nb_col);
    for (int i = 0; i<nb_lignes; i++)
    {
        for (int j = 0; j<nb_col; j++)
        {
            if (i == j)
            {
                matrice_vect.push_back(1);
            }
            else
            {
                matrice_vect.push_back(0);
            }
        }
    }
    matrice_identite.set_matrice_vect(matrice_vect);

    return matrice_identite;
}


///méthodes

vector<double> M_matrices::get_matrice_vect()
{
    return m_matrice;
}


int M_matrices::get_nb_lignes()
{
    return m_nb_lignes;
}
int M_matrices::get_nb_col()
{
    return m_nb_col;
}

void M_matrices::set_nb_lignes_col(int nb_lignes, int nb_col)
{
    m_nb_lignes = nb_lignes;
    m_nb_col = nb_col;
}

void M_matrices::set_matrice_vect(vector<double> matrice_vector)
{
    m_matrice = matrice_vector;
}

vector<int> M_matrices::donner_index_non_nul(){

    vector<int> index_non_nul;
    double matrice[m_nb_lignes][m_nb_col];

    int t = 0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice[i][j] = m_matrice[t];
            t++;
        }
    }

    for (int i=0;i<m_nb_lignes;i++){
        for(int j=0; j<m_nb_col; j++){
            if(matrice[i][j]!=0){
                index_non_nul.push_back(j+1);
                break;
            }
            else if(j+1==m_nb_col){
                index_non_nul.push_back(j+2);
            }

        }
    }

    return index_non_nul;
}

void M_matrices::swap_lignes(){

    vector<int> index_non_nul = this->donner_index_non_nul();
    double matrice[m_nb_lignes][m_nb_col];
    double new_matrice[m_nb_lignes][m_nb_col];
    double matrice_identite[m_nb_lignes][m_nb_col];
    //vector <double> matrice_identite_vect = get_matrice_identite(m_nb_lignes, m_nb_col).m_matrice;
    double new_matrice_identite[m_nb_lignes][m_nb_col];
    int t = 0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice[i][j] = m_matrice[t];
            t++;
        }
    }

    t = 0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice_identite[i][j] = m_matrice_inverse[t];
            t++;
        }
    }


    int e = 0;
    for(int i = 0; i<=m_nb_col; i++){
        for (int j = 0; j<m_nb_lignes;j++){
            if(index_non_nul[j]== i+1) {
                for (int k= 0; k<m_nb_col; k++){
                    new_matrice[e][k] = matrice[j][k];
                    if(m_nb_lignes == m_nb_col){
                        new_matrice_identite[e][k] = matrice_identite[j][k];
                    }
                }
                e++;
            }

        }
    }

    t=0;
    for(int i=0; i<m_nb_lignes; i++){
        for(int j=0; j<m_nb_col; j++){

            m_matrice[t] = new_matrice[i][j];
            m_matrice_inverse[t] = new_matrice_identite[i][j];
            t++;
        }
    }

}


void M_matrices::reduire(){

    vector<int> index_non_nul = this->donner_index_non_nul();
    double matrice_identite[m_nb_lignes][m_nb_col];
    this->swap_lignes();

    double matrice[m_nb_lignes][m_nb_col];

    int t =0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice[i][j] = m_matrice[t];
            t++;
        }
    }

    t=0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice_identite[i][j] = m_matrice_inverse[t];
            t++;
        }
    }


    double coefficient;

    for(int j = m_nb_lignes-1; j>=1;j--)
    {
        for (int l = 1; l<=j;l++){
        if (matrice[j-l][index_non_nul[j]-1] != 0 && matrice[j][index_non_nul[j]-1]!=0){

            coefficient = matrice[j-l][index_non_nul[j]-1]/matrice[j][index_non_nul[j]-1];

            for(int k = 0; k<m_nb_col; k++)
            {
                matrice[j-l][k] = matrice[j-l][k]-coefficient*matrice[j][k];

                if (m_nb_col == m_nb_lignes)
                {
                    matrice_identite[j-l][k] = matrice_identite[j-l][k]-coefficient*matrice_identite[j][k];
                }

            }
        }
        }
    }
    //}

    for (int i = 0; i<m_nb_lignes;i++){
        if(matrice[i][index_non_nul[i]-1] !=1 && matrice[i][index_non_nul[i]-1] !=0){
            double quotient = matrice[i][index_non_nul[i]-1];
            for(int k = 0; k<m_nb_col;k++){
                matrice[i][k] = matrice[i][k]/quotient;
                if(m_nb_lignes == m_nb_col){
                    matrice_identite[i][k] = matrice_identite[i][k]/quotient;
                }

            }
        }
    }

    t=0;
    for(int i=0; i<m_nb_lignes; i++){
        for(int j=0; j<m_nb_col; j++){
            m_matrice[t] = matrice[i][j];
            t++;
            }
        }
    t=0;
    for(int i=0; i<m_nb_lignes; i++){
        for(int j=0; j<m_nb_col; j++){
            m_matrice_inverse[t] = matrice_identite[i][j];
            t++;
            }
        }


}

void M_matrices::echelonner()
{

    this->swap_lignes();
    vector<int> index_non_nul = this->donner_index_non_nul();


    double matrice[m_nb_lignes][m_nb_col];
    double matrice_identite[m_nb_lignes][m_nb_col];

    int t =0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice[i][j] = m_matrice[t];
            t++;
        }
    }
    t=0;
    for (int i = 0; i< m_nb_lignes ; i++){
        for (int j = 0; j<m_nb_col; j++){
            matrice_identite[i][j] = m_matrice_inverse[t];
            t++;
        }
    }
    for(int i=0;i<m_nb_lignes-1;i++){
        index_non_nul = this->donner_index_non_nul();
        if (index_non_nul[i] == index_non_nul[i+1] && index_non_nul[i] != m_nb_col+1 ){


            double coefficient = (m_matrice[m_nb_col*(i+1)+index_non_nul[i+1]-1]) / (m_matrice[m_nb_col*(i)+index_non_nul[i]-1]);

            for(int k = 0; k<m_nb_col; k++){
            matrice[i+1][k] = matrice[i+1][k]-coefficient*matrice[i][k];
            if (m_nb_lignes == m_nb_col){
            matrice_identite[i+1][k] = matrice_identite[i+1][k]-coefficient*matrice_identite[i][k];
            }

            }

            int u=0;
            for(int i=0; i<m_nb_lignes; i++){
                for(int j=0; j<m_nb_col; j++){
                m_matrice[u] = matrice[i][j];
                m_matrice_inverse[u] = matrice_identite[i][j];
                u++;
                }
            }


        this->swap_lignes();
            u=0;
            for(int i=0; i<m_nb_lignes; i++){
                for(int j=0; j<m_nb_col; j++){
                matrice_identite[i][j] = m_matrice_inverse[u];
                matrice[i][j] = m_matrice[u];
                u++;
                }
            }

        i--;
        }
    }
}

void M_matrices::echelonner_reduire()
{
    this->echelonner();
    this->reduire();
}

void M_matrices::initialise_matrice_inverse()
{
    m_matrice_inverse = get_matrice_identite(m_nb_lignes, m_nb_col).get_matrice_vect();
}

vector<double> M_matrices::get_matrice_inverse()
{
    return m_matrice_inverse;
}

bool M_matrices::is_inversible() const
{
    M_matrices copie(*this);
    copie.echelonner_reduire();
    return (copie == get_matrice_identite(m_nb_lignes, m_nb_col));
}

M_matrices M_matrices::get_inverse()
{

    M_matrices copie_m = *this;
    copie_m.initialise_matrice_inverse();
    copie_m.echelonner();
    copie_m.reduire();
    M_matrices u(copie_m.get_nb_lignes(), copie_m.get_nb_col(), copie_m.get_matrice_inverse());

    return u;
}


//Méthodes et Opérateurs


    void M_matrices::affiche(ostream& flux) const
    {
        int t = 0;
        for (int i=0; i<m_nb_lignes; i++)
        {
            flux<<"     ||";
            for (int j=0; j<m_nb_col; j++)
            {
                flux<<" "<<m_matrice[t]<<" |";
                t++;
            }
            flux<<"|"<<endl;
        }
    }

    ostream& operator<<(ostream& flux, M_matrices const& matrice)
    {
        matrice.affiche(flux);
        return flux;
    }

    void M_matrices::recupere(){
        m_matrice.clear();
        m_matrice_inverse.clear();
        cout<<" _________________________________________"<<endl;
        cout<<"/        Entrez le nombre de lignes       \\"<<endl;
        cout<<"*******************************************"<<endl;
        cout<<"    - ";
        cin>>m_nb_lignes;
        cin.ignore();
        system("CLS");
        cout<<" _________________________________________"<<endl;
        cout<<"/       Entrez le nombre de colonnes      \\"<<endl;
        cout<<"*******************************************"<<endl;
        cout<<"    - ";
        cin>>m_nb_col;
        cin.ignore();
        system("CLS");
        this->initialise_matrice_inverse();

        double value;
        for(int i=0; i<m_nb_lignes; i++)
        {
            for (int j = 0; j<m_nb_col; j++)
            {
            cout<<" _________________________________________"<<endl;
            cout<<"/     Entrez les valeurs de la matrice    \\"<<endl;
            cout<<"*******************************************"<<endl;
            cout<<"| position ("<<i<<","<<j<<")                          |"<<endl;
            cout<<" - ";
            cin>>value;
            m_matrice.push_back(value);
            cin.ignore();
            //system("CLS");
            }
        }


    }


    istream& operator>>(istream& flux, M_matrices &matrice)// enlevé le const&
    {

        matrice.recupere();
        return flux;
    }

    //+
    M_matrices& M_matrices::operator+=(M_matrices const& autre)
    {
        try
        {
            if (autre.m_nb_lignes != m_nb_lignes || autre.m_nb_col != m_nb_col)
            {
                throw string("Operation + : impossible (taille diff)");
            }
        }
        catch(string const& str)
        {
            cerr<<str<<endl;
        }

        int t = 0;
        for (int i=0; i<m_nb_lignes; i++)
        {
            for (int j=0; j<m_nb_col; j++)
            {
                m_matrice[t] += autre.m_matrice[t];
                t++;
            }
        }

        return *this;
    }

    M_matrices operator+(M_matrices const& a, M_matrices const& b)
    {
        M_matrices copie(a);
        copie+=b;
        return copie;
    }

    //-
    M_matrices operator-(M_matrices const& a, M_matrices const& b)
    {
        M_matrices copie(a);
        copie-=b;
        return copie;
    }

    M_matrices& M_matrices::operator-=(M_matrices const& autre)
    {
        try
        {
            if (autre.m_nb_lignes != m_nb_lignes || autre.m_nb_col != m_nb_col)
            {
                throw string("Operation - : impossible (taille diff)");
            }

        }
        catch(string const& str)
        {
            cerr<<str<<endl;
        }


        int t = 0;
        for (int i=0; i<m_nb_lignes; i++)
        {
            for (int j=0; j<m_nb_col; j++)
            {
                m_matrice[t] -= autre.m_matrice[t];
                t++;
            }
        }

        return *this;
    }
    M_matrices operator*(M_matrices const& a, M_matrices const& b)
    {
        return a*b;
    }

    M_matrices operator*(M_matrices const& a, double const& lambda)
    {
        M_matrices copie(a);
        copie*=lambda;
        return copie;
    }

    M_matrices M_matrices::operator*(M_matrices const& autre)
    {
        try
        {
            if (autre.m_nb_lignes != m_nb_col)
            {
                throw string("Operation * : impossible (taille diff)");
            }
            else
            {

                M_matrices matrice_prod_M(m_nb_lignes, autre.m_nb_col);
                double matrice_prod[m_nb_lignes][autre.m_nb_col];
                double matrice_a[m_nb_lignes][m_nb_col];
                double matrice_b[autre.m_nb_lignes][autre.m_nb_col];

                int t = 0;
                for (int i = 0; i<m_nb_lignes; i++)
                {
                    for(int j = 0; j<m_nb_col; j++)
                    {
                        matrice_a[i][j] = m_matrice[t];
                        t++;
                    }
                }
                t = 0;
                for (int i = 0; i<autre.m_nb_lignes; i++)
                {
                    for(int j = 0; j<autre.m_nb_col; j++)
                    {
                        matrice_b[i][j] = autre.m_matrice[t];
                        t++;
                    }
                }
                double value = 0;
                for (int i = 0; i<m_nb_lignes; i++)
                {
                    for (int j = 0; j<autre.m_nb_col; j++)
                    {
                        for(int u = 0; u<m_nb_col; u++)
                        {
                        value += matrice_a[i][u] * matrice_b[u][j];
                        }
                        matrice_prod[i][j] = value;
                        value = 0;
                    }
                }

                matrice_prod_M.m_matrice.clear();
                for (int i = 0; i<m_nb_lignes; i++)
                {
                    for (int j = 0 ; j<autre.m_nb_col; j++)
                    {
                        matrice_prod_M.m_matrice.push_back(matrice_prod[i][j]);
                    }
                }

                return matrice_prod_M;

            }
        }
        catch(string const& str)
        {
            cerr<<str<<endl;
            return autre;
        }


    }

    M_matrices& M_matrices::operator*=(double const& lambda) //seulement pour lambda
    {
       for(int i =0; i<m_matrice.size(); i++)
       {
           m_matrice[i] *= lambda;
       }

    return *this;
    }


    M_matrices operator/(M_matrices const& a, double const& lambda) ///à définir plus tard
    {
    M_matrices copie(a);
    copie/=lambda;
    return copie;
    }

    M_matrices& M_matrices::operator/=(double const& lambda)
    {
    for(int i =0; i<m_matrice.size(); i++)
       {
           m_matrice[i] /= lambda;
       }

    return *this;
    }

    bool M_matrices::estEgal(M_matrices const& autre) const
    {
    if (m_nb_col != autre.m_nb_col || m_nb_lignes != autre.m_nb_lignes)
    {
        return false;
    }
    else {
        for (int t=0; t<m_nb_lignes * m_nb_col; t++)
        {
                if(m_matrice[t] != autre.m_matrice[t])
                {
                    return false;///voir si ça coupe le code
                }
        }
    }
        return true;
    }


    bool operator==(M_matrices const& a, M_matrices const& b)
    {
        if(a.estEgal(b))
            return true;
        else
            return false;
    }

    bool operator!=(M_matrices const& a, M_matrices const& b)
    {
        if(a.estEgal(b))
            return false;
        else
            return true;
    }


    int m_nb_col;
    int m_nb_lignes;
    vector <double> m_matrice;
    vector <double> m_matrice_inverse;

