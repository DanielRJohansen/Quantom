#include "DisplayV2.h"


DisplayV2::DisplayV2() {
#ifdef ENABLE_DISPLAY
    int success = initGLFW();
#endif
}


void DisplayV2::drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius, Int3 color) {
    float light = 0.5;
    Int3 shaded_color = color * light;
    glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);


    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); // center of circle

    for (int i = 0; i <= triangleAmount; i++) {
        light = (sin(i * 2 * PI / triangleAmount) + 1.f) / 2.f;
        shaded_color = color * light;

        glColor3ub((uint8_t)shaded_color.x, (uint8_t)shaded_color.y, (uint8_t)shaded_color.z);

        glVertex2f(
            x + (radius * cos(i * twicePi / triangleAmount)),
            y + (radius * sin(i * twicePi / triangleAmount))
        );
    }
    glEnd();
}


void DisplayV2::drawBalls(RenderBall* balls, int n_balls) {
    for (int i = 0; i < n_balls; i++) {
        RenderBall ball = balls[i];

        drawFilledCircle(ball.pos.x, ball.pos.z, ball.radius, ball.color);
    }
}

void DisplayV2::render(Simulation* simulation) {
#ifdef ENABLE_DISPLAY
    auto start = std::chrono::high_resolution_clock::now();

    RenderBall* balls = rasterizer.render(simulation);
    glClear(GL_COLOR_BUFFER_BIT);


    drawBalls(balls, rasterizer.actual_n_particles);

    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    delete[] balls;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("\tRender time: %4d ys  ", duration.count());
#endif
}

bool DisplayV2::checkWindowStatus() {
#ifdef ENABLE_DISPLAY
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }
#endif
    return true;
}

bool DisplayV2::initGLFW() {
    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(1080, 1080, "LIMA - Molecular Dynamics Engine", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    return 1;
}