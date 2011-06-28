#ifndef CORRECT_CGAL_STRUCTURE_H
#define CORRECT_CGAL_STRUCTURE_H

#include <vector>
#include <queue>


using namespace std;

class CORRECT_CGAL_STRUCTURE
{
	public:

		CORRECT_CGAL_STRUCTURE(vector< vector<float> > *_vertices_position, 
							   vector< vector<float> > *_vertices_color,
							   vector< vector<int> >   *_facets)
		{
			Vertices_position = _vertices_position;
			Vertices_color = _vertices_color;

			Facets = _facets;
		}

                vector<vector<int> > * Correct_Facet_Orientation(void)
		{

			unsigned num_vert = (unsigned)(*Vertices_position).size();
			unsigned num_faces = (unsigned)(*Facets).size();


			// 1st fuction 
			
			// to know for each vertex its corresponding face
			// numb of elements = numb of vertices;
                        vector<vector<int> > Vertices_facets;
			for (unsigned i = 0 ; i < num_vert; i++)
			{
				vector<int> Temp_vert;
				Vertices_facets.push_back(Temp_vert);
			}			
			for (unsigned i = 0; i < num_faces; i++)
			{				
				for (unsigned j = 0 ; j < (*Facets)[i].size(); j++)
				{
					int Vert_index = (*Facets)[i][j];
					Vertices_facets[Vert_index].push_back(i);
				}
			}
			
			// 2nd function
			vector<int> Status_facets;// used for 3rd function
                        vector<vector<int> > Neighboring_facets;
                        vector<vector<int> > Neighboring_edges;
			vector<int> Reoriented_facet;
			vector<int> V;
			
			for (unsigned i = 0; i < num_faces; i++)
			{				
				Neighboring_facets.push_back(V);
				Neighboring_edges.push_back(V);

				Reoriented_facet.push_back(-1);
				
				Status_facets.push_back(-1);
			}
			

			bool Check_orientation_problem = false;


			for (unsigned i = 0; i < num_faces; i++)
			{				
				for (unsigned j = 0; j < (*Facets)[i].size(); j++)
				{
					int First_vert_index = (*Facets)[i][j];
					int Second_vert_index = (*Facets)[i][(j+1) % (*Facets)[i].size()];

					for (unsigned k = 0; k < Vertices_facets[First_vert_index].size(); k++)
					{
						int Potential_neighbor_facet = Vertices_facets[First_vert_index][k];

						if ((unsigned)Potential_neighbor_facet != i )
						{
							for (unsigned l = 0; l < (*Facets)[Potential_neighbor_facet].size(); l++)
							{
								int Vertex_index = (*Facets)[Potential_neighbor_facet][l]; 
								if (Vertex_index == Second_vert_index)
								{
									Neighboring_facets[i].push_back(Potential_neighbor_facet);
									Neighboring_edges[i].push_back(First_vert_index);
									Neighboring_edges[i].push_back(Second_vert_index);									

									if (Check_orientation_problem == false)
									{
										int FV_position = -1;
										for (unsigned m = 0; m < (*Facets)[Potential_neighbor_facet].size(); m++)
										{
											if ((*Facets)[Potential_neighbor_facet][m] == First_vert_index)
											{
												FV_position = m;
												break;
											}
										}
										
										if ((FV_position - l != 1) && (l - FV_position != (*Facets)[Potential_neighbor_facet].size() - 1))
											Check_orientation_problem = true;
									}
									break;
								}
							}
						}
					}
				}				
			}			
			
			if (Check_orientation_problem == true)
			{
				// vector<int> Status_facets
				// 3rd function
				unsigned int Num_processed_facets = 0;
				std::queue<int> Queue; 
				
				int Reoriented_facet_count = 0;
				int Normal_facet_count = 0;
				do
				{
					// To find seed facet;
					for (unsigned i = 0; i < num_faces; i++)
					{
						if (Status_facets[i] == -1)
						{
							Queue.push(i);
							Reoriented_facet[i] = 0; //Good orientation
							break;
						}
					}

					while(!Queue.empty())
					{
						int Current_facet = Queue.front();
						Queue.pop();

						if (Status_facets[Current_facet] != -1)
							continue;
						
						else
						{
							Num_processed_facets++;
							Status_facets[Current_facet] = 0;
							
							for (unsigned i = 0 ; i < Neighboring_facets[Current_facet].size(); i++)
							{
								int Neighbor_index = Neighboring_facets[Current_facet][i];
								
								if ((Status_facets[Neighbor_index] == -1) && (Reoriented_facet[Neighbor_index] == -1))
								{
									int First_vert_index, Second_vert_index;

									if (Reoriented_facet[Current_facet] == 0) // current_facet is not reoriented;
									{
										First_vert_index = Neighboring_edges[Current_facet][2*i];
										Second_vert_index = Neighboring_edges[Current_facet][2*i + 1];
									}
									else // current_facet is reoriented;
									{
										First_vert_index = Neighboring_edges[Current_facet][2*i + 1];
										Second_vert_index = Neighboring_edges[Current_facet][2*i];
									}
									
									int FV_position = -1, SV_position = -1;
									
									for (unsigned j = 0; j < (*Facets)[Neighbor_index].size(); j++)
									{
										if ((*Facets)[Neighbor_index][j] == First_vert_index)
											FV_position = j;
										if ((*Facets)[Neighbor_index][j] == Second_vert_index)
											SV_position = j;																		
									}

									// Re_orient check
									// If (re_orient_needed)
									// Re_orient of neighboring facet
                                                                        //bool Is_reorientation_needed; // MT
									if ((FV_position - SV_position == 1) || (SV_position - FV_position == (int)(*Facets)[Neighbor_index].size() - 1))
									{
                                                                                //Is_reorientation_needed = false; // MT
										Reoriented_facet[Neighbor_index] = 0;

										Normal_facet_count++;
									}
									else
									{
                                                                                //Is_reorientation_needed = true; // MT
										Reoriented_facet[Neighbor_index] = 1;

										int temp = (*Facets)[Neighbor_index][0];
										(*Facets)[Neighbor_index][0] = (*Facets)[Neighbor_index][1];
										(*Facets)[Neighbor_index][1] = temp;
										Reoriented_facet_count++;		
									}								
									
									Queue.push(Neighbor_index);
								}
							}						
						}
					}

				}while(Num_processed_facets != num_faces);
			}

			return this->Facets;			
		}

	private:

		vector< vector<int> > *Facets;
		vector< vector<float> > *Vertices_position;
		vector< vector<float> > *Vertices_color;
		

};



#endif
