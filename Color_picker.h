#ifndef COLOR_PICKER_H_INCLUDED
#define COLOR_PICKER_H_INCLUDED
#include <string>
#include <iostream>
#include <stdlib.h>
#include <windows.h>


int transform_text_color(std::string color);
void set_text_color(std::string color);
void reset_text_color();
void afficher_texte_color(std::string texte, std::string color);
void set_form_color(std::string color);
void set_point_color(int color);




#endif // COLOR_PICKER_H_INCLUDED
