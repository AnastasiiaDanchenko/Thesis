#pragma once
#include "Particle.h"
#include "RigidBody.h"

#include <list>
#include <numeric>

struct GridCell2D {
	Eigen::Vector2i cellNumber;
	std::vector<Particle2D*> cellParticles;
};

struct GridCell {
	Eigen::Vector3i cellNumber;
	std::vector<Particle*> cellParticles;
};

class Grid2D {
protected:
	int gridWidth;
	int gridHeight;
	double cellSize;
	std::vector<GridCell2D> cells2D;

public:
	Grid2D();
	Grid2D(double size);

	void virtual initializeGrid();
	void virtual updateGrid();
	void neighborSearch(std::vector<Particle2D>& particles);
};

class Grid : public Grid2D {
protected:
	int gridDepth;
	std::vector<GridCell> cells;

public:
	Grid();
	Grid(double size);

	void initializeGrid() override;
	//void updateGrid() override;
	void updateGrid(std::vector<Particle*>& particles);
	void neighborSearch(std::vector<Particle*>& particles);
	void neighborSearchGhosts();
};
