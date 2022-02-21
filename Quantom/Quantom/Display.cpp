#include "Display.h"




Display::Display(Simulation* simulation) {
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    window = new sf::RenderWindow(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "LIMA DYNAMICS");    // Sometimes takes 30 sec to load, update Icue, restart pc.. I need to ditch windows... And SFML too..

    raytracer = new Raytracer(simulation, false);

    sf::Image blank_image;
    
    texture.create(WINDOW_SIZE, WINDOW_SIZE);


}



void Display::render(Simulation* simulation) {
    window->clear();
    uint8_t* image = raytracer->render(simulation);
    //texture.loadFromImage(image);

    draw(image);
    
}

void Display::animate(Trajectory* trj)
{
    uint8_t** animation = new uint8_t * [trj->n_steps];

    trj->moveToDevice();
    trj = genericMoveToDevice(trj, 1);


    //trj->n_steps = 10;  // TEMPTYEMPTEMP

    for (int i = 0; i < trj->n_steps; i++) {
        printf("Frame %d\n", i);
        animation[i] = raytracer->render(trj, i);
    }

    for (int i = 0; i < 100; i++) {
        for (int t = 0; t < trj->n_steps; t++) {
            draw(animation[t]);
            sleep(10);
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
    image = enhance(image, RAYS_PER_DIM);
    image = enhance(image, RAYS_PER_DIM*2);
    texture.update(image);
    delete[] image;
    //

    sprite.setTexture(texture, true);
    sprite.setScale(1.f, -1.f);
    sprite.setPosition(0, WINDOW_SIZE);
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


void set(uint8_t* dst, uint8_t* src) {
    for (int i = 0; i < 4; i++)
        dst[i] = src[i];
}

uint8_t* mean(uint8_t* a, uint8_t* b) {
    uint8_t* res = new uint8_t[4];
    for (int i = 0; i < 4; i++)
        res[i] = (int)((float)(a[i] + b[i]) / 2.f);
    return res;
}



uint8_t* Display::enhance(uint8_t* im, int from_size)   //doubles res
{
    uint8_t* image = new uint8_t[WINDOW_SIZE * WINDOW_SIZE * 4];
    int s = from_size;
    int b = s*2;

    for (int y = 0; y < s; y++) {
        for (int x = 0; x < s; x++) {
            set(&image[xyToIndex(x * 2, y * 2, b)], &im[xyToIndex(x, y, s)]);
        }
    }
    for (int y = 1; y < b-1; y+=2) {
        for (int x = 0; x < b-1; x+=2) {
            uint8_t* y_mean = mean(&im[xyToIndex(x/2, (y - 1)/2, s)], &im[xyToIndex(x/2, (y + 1)/2, s)]);
            set(&image[xyToIndex(x, y, b)], y_mean);
            delete[] y_mean;
        }
    }
    for (int y = 0; y < b - 1; y += 2) {
        for (int x = 1; x < b - 1; x += 2) {
            uint8_t* x_mean = mean(&im[xyToIndex((x-1) / 2, (y) / 2, s)], &im[xyToIndex((x+1) / 2, y / 2, s)]);
            set(&image[xyToIndex(x, y, b)], x_mean);
            delete[] x_mean;
        }
    }
    
    for (int y = 1; y < b - 1; y += 2) {
        for (int x = 1; x < b - 1; x += 2) {
            uint8_t* y_mean = mean(&image[xyToIndex(x, y - 1, b)], &image[xyToIndex(x, y + 1, b)]);
            uint8_t* x_mean = mean(&image[xyToIndex(x - 1, y, b)], &image[xyToIndex(x + 1, y, b)]);
            uint8_t* bi_mean = mean(y_mean, x_mean);
            set(&image[xyToIndex(x, y, b)], bi_mean);
            delete[] y_mean;
            delete[] x_mean;
            delete[] bi_mean;
        }
    }
     

    delete[] im;
    return image;
}
