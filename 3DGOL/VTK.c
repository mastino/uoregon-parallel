void writeVTK(char *** CELLS, int generation){
	FILE * fp;
	char buff[30];
	snprintf(buff, sizeof(buff), "VTK/generation_%06d.vtk", generation);
	fp = fopen(buff, "w");
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "vtk output\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET RECTILINEAR_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", ROWS+1, COLUMNS+1, DEPTH+1);
	fprintf(fp, "Z_COORDINATES %d float\n", DEPTH+1);
	for(int i=0; i <= DEPTH; i++){
		fprintf(fp, "%d ", i);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Y_COORDINATES %d float\n", COLUMNS+1);
	for(int i=0; i <= COLUMNS; i++){
		fprintf(fp, "%d ", i);
	}
	fprintf(fp, "\n");
	fprintf(fp, "X_COORDINATES %d float\n", ROWS+1);
	for(int i=0; i <= ROWS; i++){
		fprintf(fp, "%d ", i);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "CELL_DATA %d\n", ROWS*COLUMNS*DEPTH);
	fprintf(fp, "SCALARS alive int\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	int counter=0;
	for(int r=0; r < ROWS; r++){
		for(int c=0; c < COLUMNS; c++){
			for(int d=0; d < DEPTH; d++){
				counter=d+(DEPTH*(c+COLUMNS*r));
				if(counter!=0 && counter%9==0){
					fprintf(fp, "\n");
				}
				fprintf(fp, "%d ", CELLS[r][c][d]);
			}
		}
	}
}
