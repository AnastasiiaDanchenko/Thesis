#include "..\headers\Grid.h"

std::vector<Particle> particles;
std::vector<Particle2D> particles2D;
std::vector<Particle> ghostParticles;

std::vector<GridCell> grid;
std::vector<GridCell2D> grid2D;

int GRID_WIDTH;
int GRID_HEIGHT;
int GRID_DEPTH;

// Initialize uniformed grid of fluid particles
void InitFluid() {
    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
            for (int k = 0; k < PARTICLES_Z; k++) {
                Particle p;

                p.position = Eigen::Vector3d((i + 2) * SPACING, (j + 2) * SPACING, (k + 2) * SPACING);
				p.ID = particles.size();
                particles.push_back(p);
            }
        }
    }
}

void InitGhostFluid() {
    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
            Particle p;

            p.position = Eigen::Vector3d((i + 2) * SPACING, (j + 2) * SPACING, SCENE_DEPTH / 2);
            ghostParticles.push_back(p);
        }
    }
}

// Initialize boundaries
void InitBoundaries() {
    int width = WINDOW_WIDTH / SPACING - 1;
    int hight = WINDOW_HEIGHT / SPACING - 1;
    int depth = (SCENE_DEPTH / SPACING - 1) / 2;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            for (float k = 0; k < depth; k += 0.5) {
                if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1 || k < 0.5 || k > depth - 1) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * SPACING, (j + 1) * SPACING, (k + 1) * SPACING);
                    p.isFluid = false;
                    p.ID = particles.size();

                    particles.push_back(p);
                }
			}
        }
    }

    /*for (float i = 0; i < width; i ++) {
        for (float j = 0; j < hight; j ++) {
            for (float k = 0; k < depth; k ++) {
                if (i < 1 || i > width - 2 || j < 1 || j > hight - 2 || k < 1 || k > depth - 2) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * SPACING, (j + 1) * SPACING, (k + 1) * SPACING);
                    p.isFluid = false;
                    p.ID = particles.size();
                    particles.push_back(p);
                }
            }
        }
    }*/
}

void UniformGrid() {
    GRID_WIDTH = std::ceil(WINDOW_WIDTH / CELL_SIZE);
    GRID_HEIGHT = std::ceil(WINDOW_HEIGHT / CELL_SIZE);
    GRID_DEPTH = std::ceil(SCENE_DEPTH / CELL_SIZE);

    std::cout << "Using uniform grid with " << GRID_WIDTH  << "x" 
                                            << GRID_HEIGHT << "x" 
                                            << GRID_DEPTH  << " cells" << std::endl;

    grid.resize(GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH);
}

void GridUpdate() {
    // Clear grid
    for (int i = 0; i < grid.size(); i++) {
        grid[i].cellParticles.clear();
    }

    for (auto& p : particles) {
        const Eigen::Vector3i cellNumber = p.getCellNumber();
        if (cellNumber.x() < 0 || cellNumber.x() >= GRID_WIDTH || 
            cellNumber.y() < 0 || cellNumber.y() >= GRID_HEIGHT || 
            cellNumber.z() < 0 || cellNumber.z() >= GRID_DEPTH) {
            continue;
        }
        grid[cellNumber.x() + cellNumber.y() * GRID_WIDTH + 
             cellNumber.z() * GRID_WIDTH * GRID_HEIGHT].cellParticles.push_back(&p);
    }
}

void InitFluid2D() {
    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
			Particle2D p;

            /*if (j % 2 == 0) {
                p.position = Eigen::Vector2d((i + 2) * SPACING, (j + 2) * SPACING);
            }
            else {
                p.position = Eigen::Vector2d((i + 2.5) * SPACING, (j + 2) * SPACING);
            }*/

            p.position = Eigen::Vector2d((i + 2) * SPACING, (j + 2) * SPACING);
			p.ID = particles2D.size();
            particles2D.push_back(p);
		}
	}
}

void InitBoundaries2D() {
	int width = (WINDOW_WIDTH / 2) / SPACING - 1;
	int hight = (WINDOW_HEIGHT / 2) / SPACING - 1;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1) {
            //if (i < 3 || i > width - 4 || j < 3 || j > hight - 4) {
				Particle2D p;

				p.position = Eigen::Vector2d((i + 1) * SPACING, (j + 1) * SPACING);
				p.isFluid = false;
				p.ID = particles2D.size();

                particles2D.push_back(p);
			}
		}
	}
}

void InitFluidForBoundaryTest2D() {
    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
            Particle2D p;
            p.position = Eigen::Vector2d((i + 7) * SPACING, (j + 7) * SPACING);
            p.ID = particles2D.size();
            particles2D.push_back(p);
        }
    }
}

void InitMovingThroughBoundaries2D() {
    int width = (WINDOW_WIDTH / 2) / SPACING - 1;
    int hight = (WINDOW_HEIGHT / 2) / SPACING - 1;

    for (int i = 5; i < width - 5; i++) {
        for (int j = 5; j < hight - 5; j++) {
            if (i < 6 || i > width - 7 || j < 6 || j > hight - 7) {
                Particle2D p;

                p.position = Eigen::Vector2d((i + 1) * SPACING, (j + 1) * SPACING);
                p.isFluid = false;
                p.ID = particles2D.size();

                particles2D.push_back(p);
            }
        }
    }
}

void MovingBoundary2D() {
    for (int i = -((WINDOW_WIDTH / 2) / SPACING - 1) / 2; i < 0; i++) {
        for (int j = 3; j < 6; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * SPACING, j * SPACING + WINDOW_HEIGHT / 10);
            p.isFluid = false;
            p.ID = particles2D.size();
            p.velocity = Eigen::Vector2d(20, 0.0);

            particles2D.push_back(p);
        }
    }
}

Eigen::Vector2d rotatePlane(const Eigen::Vector2d& center, const Eigen::Vector2d& pos0, int angle) {
    double angleRad = angle * M_PI / 180;
    double x = center.x() + cos(angleRad) * (pos0.x() - center.x()) - sin(angleRad) * (pos0.y() - center.y());
    double y = center.y() + sin(angleRad) * (pos0.x() - center.x()) + cos(angleRad) * (pos0.y() - center.y());

    return Eigen::Vector2d(x, y);
}

void RotatingBoundary2D() {
    int width_max = WINDOW_WIDTH / SPACING / 5;
    int width_avg = WINDOW_WIDTH / SPACING / 10;

    int height_avg = WINDOW_HEIGHT / SPACING / 7;

    const Eigen::Vector2d center = Eigen::Vector2d((width_max + width_avg) * SPACING, height_avg * SPACING);
    const double angularVelocity = 0.05;

    for (int i = 0; i < width_max; i++) {
        Particle2D p0, p1, p2, p3;

        if (i == width_avg) {
            p0.position = Eigen::Vector2d((i + width_max + 2.5) * SPACING, (height_avg + 1) * SPACING);
            p0.isFluid = false;
            p0.ID = particles2D.size();

            particles2D.push_back(p0);

            p1.position = rotatePlane(center, p0.position, 45);
            p1.isFluid = false;
            p1.ID = particles2D.size();

            particles2D.push_back(p1);

            p2.position = rotatePlane(center, p0.position, 90);
            p2.isFluid = false;
            p2.ID = particles2D.size();

            particles2D.push_back(p2);

            p3.position = rotatePlane(center, p0.position, 135);
            p3.isFluid = false;
            p3.ID = particles2D.size();

            particles2D.push_back(p3);

            p0.position = Eigen::Vector2d((i + width_max - 2.5) * SPACING, (height_avg - 1) * SPACING);
        }
        else {
            p0.position = Eigen::Vector2d((i + width_max) * SPACING, height_avg * SPACING);
        }
        p0.isFluid = false;
		p0.ID = particles2D.size();

		particles2D.push_back(p0);

        p1.position = rotatePlane(center, p0.position, 45);
        p1.isFluid = false;
        p1.ID = particles2D.size();

        particles2D.push_back(p1);

        p2.position = rotatePlane(center, p0.position, 90);
        p2.isFluid = false;
        p2.ID = particles2D.size();

        particles2D.push_back(p2);

        p3.position = rotatePlane(center, p0.position, 135);
        p3.isFluid = false;
        p3.ID = particles2D.size();

        particles2D.push_back(p3);

        if (i == 0) {
            ROTATING_BOUNDARY_ID = p0.ID;
        }
	}
}

void InitFluidForRotatingTest2D() {
    for (int i = 1; i < WINDOW_WIDTH / SPACING / 4; i++) {
        for (int j = (WINDOW_HEIGHT / SPACING - 1) / 4 + 1; j < (WINDOW_HEIGHT / SPACING - 1) / 4 + 35; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * SPACING, (j + 1) * SPACING);
            p.ID = particles2D.size();

            particles2D.push_back(p);
        }
    }
}

void UniformGrid2D() {
	GRID_WIDTH = std::ceil(WINDOW_WIDTH / CELL_SIZE);
	GRID_HEIGHT = std::ceil(WINDOW_HEIGHT / CELL_SIZE);

	std::cout << "Using uniform grid with " << GRID_WIDTH << "x"
		<< GRID_HEIGHT << " cells" << std::endl;

	grid2D.resize(GRID_WIDTH * GRID_HEIGHT);
}

void GridUpdate2D() {
	// Clear grid
    for (int i = 0; i < grid2D.size(); i++) {
		grid2D[i].cellParticles.clear();
	}

    for (auto& p : particles2D) {
		const Eigen::Vector2i cellNumber = p.getCellNumber();
        if (cellNumber.x() < 0 || cellNumber.x() >= GRID_WIDTH ||
            cellNumber.y() < 0 || cellNumber.y() >= GRID_HEIGHT) {
			continue;
		}
		grid2D[cellNumber.x() + cellNumber.y() * GRID_WIDTH].cellParticles.push_back(&p);
	}
}
