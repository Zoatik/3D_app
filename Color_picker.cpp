#include <string>
#include <iostream>
#include <stdlib.h>
#include <windows.h>

using namespace std;

int transform_text_color(string color){

    int color_int = 15;//blanc

    if (color == "bleu")
    {
        color_int = 9;
    }
    else if (color == "vert")
    {
        color_int = 10;
    }
    else if (color == "rouge")
    {
        color_int = 4;
    }
    else if (color == "violet")
    {
        color_int = 13;
    }
    else if (color == "jaune")
    {
        color_int = 14;
    }
    else if (color == "blanc")
    {
        color_int = 15;
    }
    else if (color == "gris")
    {
        color_int = 8;
    }
    else
    {
        color_int = 15;
    }

    return color_int;
}

void set_text_color(string color)
{
    int color_int = transform_text_color(color);
    string sys_color = "color " + color_int;
    system(sys_color.c_str());
}

void reset_text_color()
{
   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
}

void afficher_texte_color(string texte, string color)
{
    int color_int = transform_text_color(color);
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color_int);
    cout<<texte<<endl;
    reset_text_color();
}

void set_form_color(string color)
{
    int color_int = transform_text_color(color);
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color_int);
}
void set_point_color(int color)
{
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
}

