#include "Display.h"




Display::Display(Simulation* simulation) {
    window = new sf::RenderWindow(sf::VideoMode(1000, 1000), "Quantom Simulation");


    raytracer = new Raytracer(simulation, false);

    
    sf::Image blank_image;
    blank_image.create(RAYS_PER_DIM, RAYS_PER_DIM, sf::Color(0,100, 200));
    texture.create(RAYS_PER_DIM, RAYS_PER_DIM);
    texture.loadFromImage(blank_image);




}



void Display::render(Simulation* simulation) {
    window->clear();
    uint8_t* image = raytracer->render(simulation);
    //texture.loadFromImage(image);
    texture.update(image);
    //


    delete image;
    sprite.setTexture(texture, true);

    //Flip vertically to move (0,0) from upper left corner to lower left corner.
    sprite.setScale(1.f, -1.f);
    sprite.setPosition(0,1000);

    window->draw(sprite);
    window->display();
}





void Display::terminate() {
    printf("Closing window\n");
    window->close();
}