#include "Display.h"




Display::Display() {
    window = new sf::RenderWindow(sf::VideoMode(1000, 1000), "Quantom Simulation");

    sf::CircleShape shape(100.f);
    shape.setFillColor(sf::Color::Green);


    window->clear();
    window->draw(shape);




}



void Display::render() {
    window->clear();
    window->display();
}





void Display::terminate() {
    window->close();
}