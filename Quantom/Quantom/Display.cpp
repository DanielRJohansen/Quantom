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

    draw(image);
    delete [] image;
    
}

void Display::animate(Trajectory* trj)
{
    uint8_t** animation = new uint8_t * [trj->n_steps];

    trj->moveToDevice();
    trj = genericMoveToDevice(trj, 1);


    trj->n_steps = 10;  // TEMPTYEMPTEMP

    for (int i = 0; i < trj->n_steps; i++) {
        printf("Frame %d\n", i);
        animation[i] = raytracer->render(trj, i);
    }

    for (int i = 0; i < 10; i++) {
        for (int t = 0; t < trj->n_steps; t++) {
            draw(animation[t]);
            sleep(50);
        }
    }

    printf("finished");
}





void Display::terminate() {
    printf("Closing window\n");
    window->close();
}

void Display::draw(uint8_t* image)
{
    texture.update(image);
    //


    sprite.setTexture(texture, true);

    //Flip vertically to move (0,0) from upper left corner to lower left corner.
    sprite.setScale(1.f, -1.f);
    sprite.setPosition(0, 1000);

    window->draw(sprite);
    window->display();
}

void Display::sleep(int ms)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    while (true) {
        double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
        if (duration >= ms)
            break;
    }
}
