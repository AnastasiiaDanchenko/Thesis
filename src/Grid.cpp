#include "..\headers\Grid.h"

Grid2D::Grid2D(double size) :
    gridWidth(std::ceil(parameters.windowSize.width / size)),
    gridHeight(std::ceil(parameters.windowSize.height / size)),
    cellSize(size) {
	initializeGrid();
}

void Grid2D::initializeGrid() {
	std::cout << "Using uniform grid with " << gridWidth << "x" << gridHeight << " cells" << std::endl;
    cells2D.resize(gridWidth * gridHeight);
}

void Grid2D::updateGrid() {
	for (auto& cell : cells2D) {
		cell.cellParticles.clear();
	}
    for (auto& p : particles2D) {
		const Eigen::Vector2i cellNumber = (p.position / cellSize).cast<int>();
        if (cellNumber.x() < 0 || cellNumber.x() >= gridWidth ||
            cellNumber.y() < 0 || cellNumber.y() >= gridHeight) {
            continue;
        }
        cells2D[cellNumber.x() + cellNumber.y() * gridWidth].cellParticles.push_back(&p);
    }
}

void Grid2D::neighborSearch(std::vector<Particle2D>& particles2D) {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        particles2D[i].neighbors.clear();
        const Eigen::Vector2i cellNumber = (particles2D[i].position / cellSize).cast<int>();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                if (j < 0 || j >= gridWidth || k < 0 || k >= gridHeight) { continue; }

                for (auto& p : cells2D[j + k * gridWidth].cellParticles) {
                    double distance = std::sqrt((p->position - particles2D[i].position).squaredNorm());

                    if (distance < parameters.support) {
                        particles2D[i].neighbors.push_back(p);
                    }
                }
            }
        }
    }
}

Grid::Grid(double size)
    : Grid2D(size), gridDepth(std::ceil(parameters.windowSize.depth / size)) {
    initializeGrid();
}

void Grid::initializeGrid() {
    std::cout << "Using uniform grid with " << gridWidth << "x" << gridHeight << "x" << gridDepth << " cells" 
        << std::endl;
    cells.resize(gridWidth * gridHeight * gridDepth);
}

void Grid::updateGrid() {
    for (auto& cell : cells) {
		cell.cellParticles.clear();
    }
    for (auto& p : particles) {
        const Eigen::Vector3i cellNumber = (p.position / cellSize).cast<int>();
        if (cellNumber.x() < 0 || cellNumber.x() >= gridWidth ||
            cellNumber.y() < 0 || cellNumber.y() >= gridHeight ||
            cellNumber.z() < 0 || cellNumber.z() >= gridDepth) {
			continue;
		}
		cells[cellNumber.x() + cellNumber.y() * gridWidth +
            cellNumber.z() * gridWidth * gridHeight].cellParticles.push_back(&p);
	}
}

void Grid::neighborSearch(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        particles[i].neighbors.clear();
        const Eigen::Vector3i cellNumber = (particles[i].position / cellSize).cast<int>();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                for (int l = cellNumber.z() - 1; l <= cellNumber.z() + 1; l++) {
                    if (j < 0 || j >= gridWidth || k < 0 || k >= gridHeight || l < 0 || l >= gridDepth) { continue; }

                    for (auto& p : cells[j + k * gridWidth + l * gridWidth * gridHeight].cellParticles) {
                        double distance = std::sqrt((p->position - particles[i].position).squaredNorm());

                        if (distance < parameters.support) {
                            particles[i].neighbors.push_back(p);
                        }
                    }
                }
            }
        }
    }
}

// Initialize uniformed grid of fluid particles
void InitFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            for (int k = 0; k < parameters.particlesPerDimension.z; k++) {
                Particle p;

                p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, (k + 2) * parameters.spacing);
				p.ID = particles.size();
                particles.push_back(p);
            }
        }
    }
}

void InitGhostFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle p;

            p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, parameters.windowSize.depth / 2);
            ghostParticles.push_back(p);
        }
    }
}

// Initialize boundaries
void InitBoundaries() {
    int width = parameters.windowSize.width / parameters.spacing - 1;
    int hight = parameters.windowSize.height / parameters.spacing - 1;
    int depth = (parameters.windowSize.depth / parameters.spacing - 1) / 2;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            for (float k = 0; k < depth; k += 0.5) {
                if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1 || k < 0.5 || k > depth - 1) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, (k + 1) * parameters.spacing);
                    p.isFluid = false;
                    p.ID = particles.size();
                    parameters.boundaryTestID = p.ID;

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

                    p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, (k + 1) * parameters.spacing);
                    p.isFluid = false;
                    p.ID = particles.size();
                    particles.push_back(p);
                }
            }
        }
    }*/
}

void MovingBoundary() {
    for (int i = -20; i < 0; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = 0; k < (parameters.windowSize.depth / parameters.spacing - 1) / 2; k++) {
				Particle p;

				p.position = Eigen::Vector3d(i * parameters.spacing, j * parameters.spacing, (k + 1) * parameters.spacing);
				p.isFluid = false;
				p.ID = particles.size();
				p.velocity = Eigen::Vector3d(20, 0.0, 0.0);

				particles.push_back(p);
			}
        }
    }
}

void InitFluid2D() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
			Particle2D p;

            /*if (j % 2 == 0) {
                p.position = Eigen::Vector2d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing);
            }
            else {
                p.position = Eigen::Vector2d((i + 2.5) * parameters.spacing, (j + 2) * parameters.spacing);
            }*/

            p.position = Eigen::Vector2d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing);
			p.ID = particles2D.size();
            particles2D.push_back(p);
		}
	}
}

void InitBoundaries2D() {
	int width = (parameters.windowSize.width / 2) / parameters.spacing - 1;
	int hight = (parameters.windowSize.height / 2) / parameters.spacing - 1;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1) {
            //if (i < 3 || i > width - 4 || j < 3 || j > hight - 4) {
				Particle2D p;

				p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
				p.isFluid = false;
				p.ID = particles2D.size();

                particles2D.push_back(p);
			}
		}
	}
}

void InitFluidForBoundaryTest2D() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle2D p;
            p.position = Eigen::Vector2d((i + 7) * parameters.spacing, (j + 7) * parameters.spacing);
            p.ID = particles2D.size();
            particles2D.push_back(p);
        }
    }
}

void InitMovingThroughBoundaries2D() {
    int width = (parameters.windowSize.width / 2) / parameters.spacing - 1;
    int hight = (parameters.windowSize.height / 2) / parameters.spacing - 1;

    for (int i = 5; i < width - 5; i++) {
        for (int j = 5; j < hight - 5; j++) {
            if (i < 6 || i > width - 7 || j < 6 || j > hight - 7) {
                Particle2D p;

                p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
                p.isFluid = false;
                p.ID = particles2D.size();

                particles2D.push_back(p);
            }
        }
    }
}

void MovingBoundary2D() {
    for (int i = -((parameters.windowSize.width / 2) / parameters.spacing - 1) / 2; i < 0; i++) {
        for (int j = 3; j < 6; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * parameters.spacing, j * parameters.spacing + parameters.windowSize.height / 10);
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
    int width_max = parameters.windowSize.width / parameters.spacing / 5;
    int width_avg = parameters.windowSize.width / parameters.spacing / 10;

    int height_avg = parameters.windowSize.height / parameters.spacing / 7;

    const Eigen::Vector2d center = Eigen::Vector2d((width_max + width_avg) * parameters.spacing, height_avg * parameters.spacing);
    const double angularVelocity = 0.05;

    for (int i = 0; i < width_max; i++) {
        Particle2D p0, p1, p2, p3;

        if (i == width_avg) {
            p0.position = Eigen::Vector2d((i + width_max + 2.5) * parameters.spacing, (height_avg + 1) * parameters.spacing);
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

            p0.position = Eigen::Vector2d((i + width_max - 2.5) * parameters.spacing, (height_avg - 1) * parameters.spacing);
        }
        else {
            p0.position = Eigen::Vector2d((i + width_max) * parameters.spacing, height_avg * parameters.spacing);
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
            parameters.boundaryTestID = p0.ID;
        }
	}
}

void InitFluidForRotatingTest2D() {
    for (int i = 1; i < parameters.windowSize.width / parameters.spacing / 4; i++) {
        for (int j = (parameters.windowSize.height / parameters.spacing - 1) / 4 + 1; j < (parameters.windowSize.height / parameters.spacing - 1) / 4 + 35; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
            p.ID = particles2D.size();

            particles2D.push_back(p);
        }
    }
}
