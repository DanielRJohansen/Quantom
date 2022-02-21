#include "DisplayV2.h"


DisplayV2::DisplayV2() {
    bool success = initGLFW();
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

void drawTwoSquares() {
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3f(1.0, 0.25, 1.0);
    glBegin(GL_POLYGON);
    glVertex3f(-0.5, -0.50, 0.0);
    glVertex3f(0.75, 0.25, .50);
    glVertex3f(0.75, 0.75, .50);
    glVertex3f(0.25, 0.75, .50);
    glEnd();


    glColor3f(.0, 0.25, 1.0);
    glBegin(GL_POLYGON);
    glVertex3f(0.1, 0.1, 0.40);
    glVertex3f(0.5, 0.1, 0.40);
    glVertex3f(0.5, 0.5, 0.40);
    glVertex3f(0.1, 0.5, 0.40);
    glEnd();
}

void DisplayV2::drawBalls(RenderBall* balls, int n_balls) {
    for (int i = 0; i < n_balls; i++) {
        RenderBall ball = balls[i];

        drawFilledCircle(ball.pos.x, ball.pos.z, ball.radius, ball.color);
    }
}

void DisplayV2::render(Simulation* simulation) {
    auto start = std::chrono::high_resolution_clock::now();

    RenderBall* balls = rasterizer.render(simulation);
    glClear(GL_COLOR_BUFFER_BIT);


    drawBalls(balls, rasterizer.actual_n_particles);

    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("\tRender time: %4d ys  ", duration.count());

}

bool DisplayV2::checkWindowStatus() {
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }
        
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
}