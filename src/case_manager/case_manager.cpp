#include <case_manager.hpp>
 
using namespace std;

// unions for supported settings

typedef std::vector<ryml::csubstr> rymlKeys;



// helper functions from ryml library

C4_SUPPRESS_WARNING_MSVC_WITH_PUSH(4996) // fopen: this function or variable may be unsafe
/** load a file from disk and return a newly created CharContainer */
template<class CharContainer>
size_t file_get_contents(const char *filename, CharContainer *v)
{
   ::FILE *fp = ::fopen(filename, "rb");
   C4_CHECK_MSG(fp != nullptr, "could not open file");
   ::fseek(fp, 0, SEEK_END);
   long sz = ::ftell(fp);
   v->resize(static_cast<typename CharContainer::size_type>(sz));
   if(sz)
   {
       ::rewind(fp);
       size_t ret = ::fread(&(*v)[0], 1, v->size(), fp);
       C4_CHECK(ret == (size_t)sz);
   }
   ::fclose(fp);
   return v->size();
}

/** load a file from disk into an existing CharContainer */
template<class CharContainer>
CharContainer file_get_contents(const char *filename)
{
   CharContainer cc;
   file_get_contents(filename, &cc);
   return cc;
}


CaseManager::CaseManager(std::string& caseFileName, VisItDataCollection& rdc)
:
   contents(""),
   restart(false),
   restartCycle(0),
   nSavedFrames(0),
   restart_data_c (rdc),
   restart_current_cycle_name(""),
   restart_deleted_cycle_name(""),
   serRefLevels(0),
   parRefLevels(0),
   adaptive_mesh(false),
   spatialOrder(1),
   total_error_fraction(0.0),
   max_elem_error(0.0),
   hysteresis(0.0),
   nc_limit(0),
   prefer_conforming_refinement(false),
   icType("constant"),
   rsolverType("Rusanov"),
   ICInterface(NULL),
   origin(3),
   normal(3),
   c1min(3),
   c1max(3),
   c2min(3),
   c2max(3),
   a(NULL),
   b(NULL),
   c(NULL),
   checkTotalEnergy(false),
   writeIndicators(false),
   linearize(false),
   haveLastHope(true)

{
   caseDir = std::filesystem::current_path().string();
   caseFileName = caseDir + "/" + caseFileName;
   std::string contents = file_get_contents<std::string>(caseFileName.c_str());
   settings = new ryml::Tree(ryml::parse_in_arena(ryml::to_csubstr(contents))); 
}

CaseManager::~CaseManager()
{
   delete[] c;
   delete[] b;
   delete[] a;
   delete ICInterface;
   delete settings;
   // delete contents;
   origin.SetSize(1);
   origin[0] = 0;
   normal.SetSize(1);
   normal[0] = 0;
};

void CaseManager::loadRestartSettings()
{
   restart = readOrDefault<bool>(settings,rymlKeys({"time","restart"}), 0);
   if (myRank == 0) std::cout << "Case is restarted: " << restart << std::endl;

   restartCycle = read<int>(settings,rymlKeys({"time","restartCycle"}));
   nSavedFrames = read<int>(settings,rymlKeys({"time","nSavedFrames"}));
   
}

void CaseManager::loadPhysics()
{
   ryml::csubstr phType = read<ryml::csubstr>(settings,rymlKeys({"physics","type"}));
   if (myRank == 0) std::cout << "Gas EoS type: " << phType << std::endl;

   specific_heat_ratio = readOrDefault<double>(settings,rymlKeys({"physics","gamma"}),1.4);
   if (myRank == 0) std::cout << "* Specific heat ratio: " << specific_heat_ratio << std::endl;
   
   gas_constant = readOrDefault<double>(settings,rymlKeys({"physics","R"}),1.0);
   if (myRank == 0) std::cout << "* Gas constant: " << gas_constant << std::endl;
   
   if (phType == "covolume")
   {
       covolume_constant = read<double>(settings,rymlKeys({"physics","covolume"})   );
       if (myRank == 0) std::cout << "* Covolume constant: " << covolume_constant << std::endl;
   }
   else
       covolume_constant = 0.0;
}

void CaseManager::loadSpatialAccuracyOrder()
{
   spatialOrder = read<int>(settings,rymlKeys({"spatial","polyOrder"}));
   if (myRank == 0) std::cout << "Order of polynoms: " << spatialOrder << std::endl;
}

void CaseManager::loadPostprocessSettings()
{
   checkTotalEnergy = readOrDefault<bool>(settings,rymlKeys({"postProcess","checkTotalEnergy"}), 0);
   if (myRank == 0) std::cout << "Check total energy: " << checkTotalEnergy << std::endl;

   writeIndicators = readOrDefault<bool>(settings,rymlKeys({"postProcess","writeIndicators"}), 0);
   if (myRank == 0) std::cout << "Write indicators field: " << writeIndicators << std::endl;
}

void CaseManager::loadInitialConditions()
{
   icType = read<ryml::csubstr>(settings,rymlKeys({"internalField","type"}));
   if (myRank == 0) std::cout << "Type of problem: " << icType << std::endl;

   if (icType == "constant")
   {
       Vector sol(num_equation);

       sol[0] = read<double>(settings,rymlKeys({"internalField","rho"}));
       sol[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),0);
       sol[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),1);
       sol[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),2);
       sol[num_equation-1] = read<double>(settings,rymlKeys({"internalField","p"}));

       ICInterface = new ICConstant(sol);
   }
   else if (icType == "planeBreakup")
   {
       Vector sol1(num_equation);
       Vector sol2(num_equation);

       origin[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),0);
       origin[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),1);
       origin[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),2);

       normal[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","normal"}),0);
       normal[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","normal"}),1);
       normal[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","normal"}),2);

       sol1[0] = read<double>(settings,rymlKeys({"internalField","left","rho"}));
       sol1[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),0);
       sol1[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),1);
       sol1[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),2);
       sol1[num_equation-1] = read<double>(settings,rymlKeys({"internalField","left","p"}));
       
       sol2[0] = read<double>(settings,rymlKeys({"internalField","right","rho"}));
       sol2[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),0);
       sol2[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),1);
       sol2[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),2);
       sol2[num_equation-1] = read<double>(settings,rymlKeys({"internalField","right","p"}));
       
       ICInterface = new ICPlaneBreakup(sol1, sol2, origin, normal);
   }
   else if (icType == "sphericalBreakup")
   {
       Vector sol1(num_equation);
       Vector sol2(num_equation);
       double radius = 0.0;

       origin[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),0);
       origin[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),1);
       origin[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),2);

       radius = read<double>(settings,rymlKeys({"internalField","radius"}));

       sol1[0] = read<double>(settings,rymlKeys({"internalField","left","rho"}));
       sol1[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),0);
       sol1[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),1);
       sol1[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),2);
       sol1[num_equation-1] = read<double>(settings,rymlKeys({"internalField","left","p"}));
       
       sol2[0] = read<double>(settings,rymlKeys({"internalField","right","rho"}));
       sol2[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),0);
       sol2[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),1);
       sol2[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),2);
       sol2[num_equation-1] = read<double>(settings,rymlKeys({"internalField","right","p"}));
       
       ICInterface = new ICSphericalBreakup(sol1, sol2, origin, radius);
   }
   else if (icType == "densityPulse")
   {
       Vector sol(num_equation);
       double epsilon = 0.0;
       
       origin[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),0);
       origin[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),1);
       origin[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","origin"}),2);

       epsilon = read<double>(settings,rymlKeys({"internalField","epsilon"}));

       sol[0] = 0.0;
       sol[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),0);
       sol[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),1);
       sol[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","U"}),2);
       sol[num_equation-1] = read<double>(settings,rymlKeys({"internalField","p"}));

       ICInterface = new ICDensityPulse(sol, origin, epsilon);
   }
   else if (icType == "twoHexahedronsBreakup")
   {
       Vector sol1(num_equation);
       Vector sol2(num_equation);

       c1min[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","min1"}),0);
       c1min[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","min1"}),1);
       c1min[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","min1"}),2);

       c1max[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","max1"}),0);
       c1max[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","max1"}),1);
       c1max[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","max1"}),2);

       c2min[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","min2"}),0);
       c2min[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","min2"}),1);
       c2min[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","min2"}),2);

       c2max[0] = readVectorComponent<double>(settings,rymlKeys({"internalField","max2"}),0);
       c2max[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","max2"}),1);
       c2max[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","max2"}),2);

       sol1[0] = read<double>(settings,rymlKeys({"internalField","left","rho"}));
       sol1[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),0);
       sol1[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),1);
       sol1[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","left","U"}),2);
       sol1[num_equation-1] = read<double>(settings,rymlKeys({"internalField","left","p"}));
       
       sol2[0] = read<double>(settings,rymlKeys({"internalField","right","rho"}));
       sol2[1] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),0);
       sol2[2] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),1);
       sol2[3] = readVectorComponent<double>(settings,rymlKeys({"internalField","right","U"}),2);
       sol2[num_equation-1] = read<double>(settings,rymlKeys({"internalField","right","p"}));
       
       ICInterface = new ICTwoHexahedronsBreakup(sol1, sol2, c1min, c1max, c2min, c2max);
   }

   else
   {
       if (myRank == 0)
       {
            std::cout << "Wrong type of initial condition" << std::endl
              << "Available types: "
              << "constant, planeBreakup, sphericalBreakup, twoHexahedronsBreakup" << std::endl;
       }
   }

   if (myRank == 0) std::cout << "Initial conditions OK" << std::endl;
}

void CaseManager::loadMesh(ParMesh*& pmesh)
{
   // resh mesh filename and convert it to char* 

   ryml::csubstr meshFileRaw = read<ryml::csubstr>(settings,rymlKeys({"mesh","file"}));
   std:string meshFileStr(meshFileRaw.str, meshFileRaw.len);
   meshFile = meshFileStr.c_str();
   
   // read mesh settings

   serRefLevels = readOrDefault<int>(settings,rymlKeys({"mesh","serialRefLevels"}), 0);
   parRefLevels = readOrDefault<int>(settings,rymlKeys({"mesh","parallelRefLevels"}), 0);
   adaptive_mesh = readOrDefault<bool>(settings,rymlKeys({"mesh","adaptive"}), 0);
   local_refinement = readOrDefault<bool>(settings,rymlKeys({"mesh","localRefinement"}), 0);

   if (myRank == 0)
   {
       std::cout << "Reading mesh from " << meshFile << "..." << std::endl;
       std::cout << "* Serial refinement levels = " << serRefLevels << std::endl;
       std::cout << "* Parallel refinement levels = " << parRefLevels << std::endl;
       std::cout << "* Adaptive: " << (adaptive_mesh ? "true" : "false") << std::endl;
       std::cout << "* Local refinement: " << (local_refinement ? "true" : "false") << std::endl;
   }
   // read adaptive mesh settings
//?//

   if (restart) // load par mesh in case of restart
   {
       // 2.1. Read the given data collection.
       restart_data_c.Load(restartCycle);
       if (restart_data_c.Error())
       {
            if (myRank == 0)
            {
                std::cout << "Error loading data collection" << std::endl;
            }
            return;
       }
       // 2.2. Get par mesh pieces from the loaded collection
       pmesh = dynamic_cast<ParMesh*>(restart_data_c.GetMesh());
       if (!pmesh)
       {
            if (myRank == 0)
            {
                std::cout << "The given data collection does not have a parallel mesh."
                        << std::endl;
            }
            return;
       }
   }
   else //load mesh from file and make a partition
   {
       // 2.1. Read the serial mesh on all processors, refine it in serial, then
       //    partition it across all processors and refine it in parallel.

       std::ifstream meshStream;
       meshStream.open(meshFile);

       Mesh *mesh = new Mesh(meshStream, 1, 1);
       mesh->EnsureNCMesh(true);

       for (int l = 0; l < serRefLevels; l++)
       {
            mesh->UniformRefinement();
       }

       if (local_refinement)
       {
            localRefLevels = readOrDefault<int>(settings,rymlKeys({"mesh","localRefinement","refLevels"}), 0);
            localRefInside = readOrDefault<bool>(settings,rymlKeys({"mesh","localRefinement","inside"}), 0);
            localRefType = read<ryml::csubstr>(settings,rymlKeys({"mesh","localRefinement","type"}));

            if (localRefType == "sphere")
            {
                double radius = 0.0;
                
                origin[0] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","origin"}),0);
                origin[1] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","origin"}),1);
                origin[2] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","origin"}),2);

                radius = read<double>(settings,rymlKeys({"mesh","localRefinement","radius"}));
                
                SphericalLocalRefiner locref(origin,radius,mesh,localRefLevels,localRefInside);
                locref.Refine();
            }
            else if (localRefType == "hex")
            {
                c1min[0] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","min1"}),0);
                c1min[1] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","min1"}),1);
                c1min[2] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","min1"}),2);

                c1max[0] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","max1"}),0);
                c1max[1] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","max1"}),1);
                c1max[2] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","max1"}),2);

                HexahedronLocalRefiner locref(c1min,c1max,mesh,localRefLevels,localRefInside);
                locref.Refine();
            }
            else if (localRefType == "cylinder")
            {
                double radius = 0.0;
                
                c1min[0] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p1"}),0);
                c1min[1] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p1"}),1);
                c1min[2] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p1"}),2);

                c1max[0] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p2"}),0);
                c1max[1] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p2"}),1);
                c1max[2] = readVectorComponent<double>(settings,rymlKeys({"mesh","localRefinement","p2"}),2);
                
                radius = read<double>(settings,rymlKeys({"mesh","localRefinement","radius"}));
                
                CylindricalLocalRefiner locref(c1min,c1max,radius,mesh,localRefLevels,localRefInside);
                locref.Refine();
            }
            else
            {
                if (myRank == 0) 
                     std::cout << "Unknown local refiner type: " << localRefType << '\n';
                exit(1);
            }
       }



       if (myRank == 0) { std::cout << "Number of cells: " << mesh->GetNE() << std::endl; }

       // 2.2. Define a parallel mesh by a partitioning of the serial mesh. Refine
       //    this mesh further in parallel to increase the resolution. Once the
       //    parallel mesh is defined, the serial mesh can be deleted.
       pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
       delete mesh;

       for (int l = 0; l < parRefLevels; l++)
       {
            pmesh->UniformRefinement();
       }
       restart_data_c.SetFormat(1);
       restart_data_c.SetMesh(pmesh);  
   } 

   if (myRank == 0) std::cout << "Dimension = " << pmesh->Dimension() << std::endl;

   spaceDim = pmesh->SpaceDimension();
   readPhysicalNames(pmesh);
   // minimizeAttributes(pmesh);
}

void CaseManager::readPhysicalNames(ParMesh*& mesh)
{
   string buff;
   int number_of_physical_groups;
   int dim_r;
   std::string phys_group_name_r;
   int tag_r;

   std::map<std::string,int>& map_phys_names_tag = (mesh->SpaceDimension() == 2) ? map_phys_names_tag_2D : map_phys_names_tag_3D;
   std::map<std::string,int>& map_bdr_names_tag  = (mesh->SpaceDimension() == 2) ? map_phys_names_tag_1D : map_phys_names_tag_2D;

   std::ifstream input(meshFile);

   while (input >> buff)
   {
       if (buff == "$PhysicalNames") 
       {
            input >> number_of_physical_groups;
            getline(input, buff);

            for (int iGroup = 0; iGroup < number_of_physical_groups; ++iGroup)
            {
                input >> dim_r >> tag_r >> phys_group_name_r;

                if ( phys_group_name_r.front() == '"' ) 
                {
                     phys_group_name_r.erase( 0, 1 ); // erase the first character
                     phys_group_name_r.erase( phys_group_name_r.size() - 1 ); // erase the last character
                }

                std::cout << phys_group_name_r << std::endl;

                if ( (dim_r == 1 && mesh->SpaceDimension() == 2) || (dim_r == 2 && mesh->SpaceDimension() == 3) )
                {
                     map_bdr_names_tag[phys_group_name_r] = tag_r;
                }
                else if ( (dim_r == 2 && mesh->SpaceDimension() == 2) || (dim_r == 3 && mesh->SpaceDimension() == 3) )
                {
                     map_phys_names_tag[phys_group_name_r] = tag_r;
                }
                else
                {
                     std::cout << "Wrong matching of dimension for physical group"
                          << "(" << phys_group_name_r << ":" << dim_r << ") "
                          << "and space dimension " 
                          << mesh->SpaceDimension()
                          << std::endl;
                     exit(1);
                }
            }

            break;
       }
       else
       if (buff == "$EndElements")
       {
            std::cout << "No physical boundaries in this mesh";
            break;
       }
   }

   input.close();
}

void CaseManager::minimizeAttributes(ParMesh*& mesh)
{
   // start from internal cells; default attribute should be 1
   mesh->bdr_attributes.Print(std::cout << "mesh.bdr_attributes = ");
   mesh->attributes.Print(std::cout << "mesh.attributes = ");

   // edit attributes for every cell
   for (int iCell = 0; iCell < mesh->GetNE(); ++iCell)
   {
       int curAttr = mesh->GetAttribute(iCell);
       for (int iAttr = 0; iAttr < mesh->attributes.Size(); ++iAttr)
       {
            if (curAttr == mesh->attributes[iAttr])
            {
                mesh->SetAttribute(iCell, iAttr + 1);
                break;
            }
       }
   }

   // edit attributes inside map for names
   std::map<std::string,int>& map_phys_names_tag = (mesh->SpaceDimension() == 2) ? map_phys_names_tag_2D : map_phys_names_tag_3D;

   for(std::map<std::string,int>::iterator it = map_phys_names_tag.begin(); it != map_phys_names_tag.end(); it++)
   {
       for (int iAttr = 0; iAttr < mesh->attributes.Size(); ++iAttr)
       {
            if (it->second == mesh->attributes[iAttr])
            {
                it->second = iAttr + 1;
                break;
            }
       }
   }

   // edit mesh attributes
   for (int i = 0; i < mesh->attributes.Size(); ++i)
       mesh->attributes[i] = i+1;

   // continue with boundary attributes, start attr numeration after cells attributer=s
   bool hasUnitBdrAttribute = 0;

   //edit attribute for every boundary element
   for (int iCell = 0; iCell < mesh->GetNBE(); ++iCell)
   {
       int curAttr = mesh->GetBdrAttribute(iCell);
       for (int iAttr = 0; iAttr < mesh->bdr_attributes.Size(); ++iAttr)
       {
            // if (curAttr == 1)
            // {
            //    hasUnitBdrAttribute = 1;
            //    break;
            // }
            if (curAttr == mesh->bdr_attributes[iAttr])
            {
                mesh->SetBdrAttribute(iCell, iAttr + 1);
                break;
            }
       }
   }

   // std::cout << "hasUnitBdrAttribute = " << hasUnitBdrAttribute << std::endl;

   // edit boundary attributes inside map for names
   std::map<std::string,int>& map_bdr_names_tag = (mesh->SpaceDimension() == 2) ? map_phys_names_tag_1D : map_phys_names_tag_2D;

   for(std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
   {
       for (int iAttr = 0; iAttr < mesh->bdr_attributes.Size(); ++iAttr)
       {
            if (it->second == mesh->bdr_attributes[iAttr])
            {
                it->second = iAttr + 1;
                break;
            }
       }
   }

   // edit mesh attributes
   // for (int i = 0; i < mesh->bdr_attributes.Size() - hasUnitBdrAttribute; ++i)
   //    mesh->bdr_attributes[i + hasUnitBdrAttribute] = i + 1 + mesh->attributes.Size();
   
   for (int i = 0; i < mesh->bdr_attributes.Size(); ++i)
   {
       mesh->bdr_attributes[i] = i + 1;
   }

   // if (hasUnitBdrAttribute)
   //    mesh->bdr_attributes[0] = 1;

}

void CaseManager::loadAdaptiveMeshSettings(ThresholdRefiner& refiner,ThresholdDerefiner& derefiner)
{
   // c4::from_chars(settings["mesh"]["totalErrorFraction"].val(), &total_error_fraction);
   // c4::from_chars(settings["mesh"]["maxElemError"].val(), &max_elem_error);
   // c4::from_chars(settings["mesh"]["hysteresis"].val(), &hysteresis);
   // c4::from_chars(settings["mesh"]["preferConformingRefinement"].val(), &prefer_conforming_refinement);
   // c4::from_chars(settings["mesh"]["nonConformingLimit"].val(), &nc_limit);

   // refiner.SetTotalErrorFraction(total_error_fraction); // use purely local threshold

   // refiner.SetLocalErrorGoal(max_elem_error);
   // if (prefer_conforming_refinement)
   //    refiner.PreferConformingRefinement();
   // refiner.SetNCLimit(nc_limit);

   // derefiner.SetThreshold(hysteresis * max_elem_error);
   // derefiner.SetNCLimit(nc_limit);
}

void CaseManager::checkNumEqn(int dim)
{
   if (dim == 2)
   {
       num_equation = 4;
   }
   else if (dim == 3)
   {
       num_equation = 5;
   }
   else
   {
       if (myRank == 0)
       {
            std::cout << "Wrong dimension " << dim << std::endl;
            std::cout << "Please check your mesh file" << std::endl;
       }
       exit(1);
   }

   if (myRank == 0) std::cout << "Number of equations = " << num_equation << std::endl;
}

void CaseManager::loadInitialSolution( ParFiniteElementSpace& vfes, const Array<int>& offsets, BlockVector& u_block, ParGridFunction& sol)
{
   ParGridFunction *saved_sol = NULL;

   // In case of restart read solution from the VisIt collection
   if (restart)
   {
       saved_sol = restart_data_c.GetParField("restart_conservative_vars");
       if (!saved_sol)
       {
            if (myRank == 0)
            {
                std::cout << "The given data collection has no 'restart_conservative_vars' field."
                        << std::endl;
                exit(2);
            }
       }
       else
       {
            u_block.Update(*saved_sol,offsets);
            sol.MakeRef(&vfes, u_block, 0);
       }
   }
   else
   {
       // Initialize the state.
       VectorFunctionCoefficient u0(num_equation, ICInterface->setIC());
       sol.ProjectCoefficient(u0);
       
       // Save the state to the VisIt
       restart_data_c.RegisterField("restart_conservative_vars",&sol);
       restart_data_c.Save();
   }
}

void CaseManager::loadRiemannSolver(RiemannSolver*& rsolver)
{
   //read riemann solver
   rsolverType = read<ryml::csubstr>(settings,rymlKeys({"spatial","riemannSolver"}));
   if (myRank == 0) std::cout << "Riemann solver type: " << rsolverType << std::endl;

   if (rsolverType == "Rusanov")
   {
       rsolver = new RiemannSolverRusanov();
   }
   else if (rsolverType == "LLF")
   {
       rsolver = new RiemannSolverLLF();
   }
   else if (rsolverType == "HLL")
   {
       rsolver = new RiemannSolverHLL();
   }
   else if (rsolverType == "HLLC")
   {
       rsolver = new RiemannSolverHLLC();
   }
   else 
   {
       if (myRank == 0) std::cout << "Unknown Riemann solver type: " << rsolverType << '\n';
       exit(1);
   }
}

int CaseManager::getFinDiffGroupAttribute(ryml::csubstr fdGroupName)
{
   std::map<std::string,int>& map_phys_names_tag = (spaceDim == 2) ? map_phys_names_tag_2D : map_phys_names_tag_3D;

   for(std::map<std::string,int>::iterator it = map_phys_names_tag.begin(); it != map_phys_names_tag.end(); it++)
   {
       if (it->first == fdGroupName)
       {
            return it->second;
       }
   }

   return -1;
}

void CaseManager::loadLimiter(Averager& avgr, Indicator*& ind, Limiter*& l, const Array<int>& offsets, const int dim, ParFiniteElementSpace& fes_const, ParFiniteElementSpace& vfes)
{
   ryml::csubstr limiterType = read<ryml::csubstr>(settings,rymlKeys({"spatial","limiter","type"}));
   if (myRank == 0) std::cout << limiterType << std::endl;

   linearize = readOrDefault<bool>(settings,rymlKeys({"spatial","limiter","linearize"}), 0);
   haveLastHope = readOrDefault<bool>(settings,rymlKeys({"spatial","limiter","haveLastHope"}), 1);
   if (myRank == 0) std::cout << "Additional linearization: " << linearize << std::endl;
   if (myRank == 0) std::cout << "Last hope limiter (cut slopes in cell in case of nonphysical values in vertices after limiting): " << haveLastHope << std::endl;
   
   ryml::csubstr fdGroupName = readOrDefault<ryml::csubstr>(settings,rymlKeys({"spatial","limiter","cutSlopeGroup"}), "");

   if (fdGroupName != "")
   {
       fdGroupAttr = getFinDiffGroupAttribute(fdGroupName);
       if (myRank == 0) std::cout << "Cut slopes in the following group of cells: " << fdGroupName << std::endl;
   }
   else
   {
       fdGroupAttr = -1;
   }
   
   ryml::csubstr indicatorType = read<ryml::csubstr>(settings,rymlKeys({"spatial","limiter","indicator"}));

   // Define indicator object

   if (indicatorType == "Nowhere")
   {
       ind = new IndicatorNowhere(avgr, &vfes, &fes_const, offsets, dim);
   }
   else if (indicatorType == "Everywhere")
   {
       ind = new IndicatorEverywhere(avgr, &vfes, &fes_const, offsets, dim); 
   }
   else if (indicatorType == "BJ")
   {
       ind = new IndicatorBJ(avgr, &vfes, &fes_const, offsets, dim); 
   }
   else if (indicatorType == "Venkatakrishnan" || indicatorType == "Venkata")
   {
       ind = new IndicatorVenkatakrishnan(avgr, &vfes, &fes_const, offsets, dim); 
   }
   else if (indicatorType == "Michalak")
   {
       ind = new IndicatorMichalak(avgr, &vfes, &fes_const, offsets, dim); 
   }
   else if (indicatorType == "Shu")
   {
       ind = new IndicatorShu(avgr, &vfes, &fes_const, offsets, dim);   
   }
   else
   {
       if (myRank == 0) std::cout << "Unknown indicator type: " << indicatorType << '\n';
       exit(1);
   }  

   // Define limiter object
   
   if (limiterType == "FinDiff")
   {
       l = new LimiterFinDiff(*ind, avgr, &vfes, offsets, linearize, haveLastHope, fdGroupAttr, dim);
   }
   else if (limiterType == "Multi")
   {
       l = new LimiterMultiplier(*ind, avgr, &vfes, offsets, linearize, haveLastHope, fdGroupAttr, dim);
   }
   else if (limiterType == "None")
   {
       l = new LimiterNone(*ind, avgr, &vfes, offsets, linearize, haveLastHope, fdGroupAttr, dim);
   }
   else
   {
       if (myRank == 0) std::cout << "Unknown limiter type: " << limiterType << '\n';
       exit(1);
   }

   if (myRank == 0) std::cout << "Indicator type: " << indicatorType << "; Limiter type: " << limiterType << std::endl;
}

void CaseManager::addBoundaryIntegrators(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim)
{
   int tag = 0;
   ryml::csubstr nameYaml;
   ryml::csubstr type;
    
   // for fully periodic meshes
   if (mesh.bdr_attributes.Size() == 0) 
       return;

   std::map<std::string,int>& map_phys_names_tag = (mesh.SpaceDimension() == 2) ? map_phys_names_tag_2D : map_phys_names_tag_3D;
   std::map<std::string,int>& map_bdr_names_tag  = (mesh.SpaceDimension() == 2) ? map_phys_names_tag_1D : map_phys_names_tag_2D;

   if (myRank == 0) 
   {
       mesh.bdr_attributes.Print(std::cout << "mesh.bdr_attributes = ");
       mesh.attributes.Print(std::cout << "mesh.attributes = ");

       std::cout << "Boundary mapping:\n";
       for ( std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
       {
            std::cout << it->first << ' ' << it->second << std::endl;
       }
       std::cout << "Physical group mapping:\n";
       for ( std::map<std::string,int>::iterator it = map_phys_names_tag.begin(); it != map_phys_names_tag.end(); it++)
       {
            std::cout << it->first << ' ' << it->second << std::endl;
       } 
   }


   bdr_markers.resize(mesh.bdr_attributes.Max());
   for (int i = 0; i < mesh.bdr_attributes.Max(); ++i)
       bdr_markers[i].SetSize(mesh.bdr_attributes.Max());

   if (myRank == 0) mesh.bdr_attributes.Print(std::cout << "Boundary attributes in mesh:\n");   

   
   if ( ((*settings)["boundaryField"].num_children() != mesh.bdr_attributes.Size()) )
   {
       if (myRank == 0)
       {
            std::cout << "=====\nPossible mismatch of boundary conditions and geometric boundary groups." << std::endl;
            std::cout << "settings[boundaryField].num_children() = " << (*settings)["boundaryField"].num_children() << std::endl;
            std::cout << "mesh.bdr_attributes.Size() = " << mesh.bdr_attributes.Size() << std::endl;
            std::cout << "Boundary names:" << std::endl;
            for (ryml::NodeRef node : (*settings)["boundaryField"].children())
            {
                std::cout << " * " << node.key() << std::endl;
            }
            std::cout << "Geometric boundaries:\n";
            for ( std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
            {
                std::cout << " * " << it->first << std::endl;
            }
       }
       
       exit(1);
   }

   if (myRank == 0) std::cout << "Boundary conditions:" << std::endl;

   // read patch names and tags from mesh file $PhysicalNames
   int patchNumber = 0;
   for (auto node : (*settings)["boundaryField"].children())
   {
       // find tag for the name
       tag = -1;
       std::cout << node.key() << std::endl;
       nameYaml = node.key();
       std::string name (nameYaml.str, nameYaml.len);
       for ( std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
       {
            if (name == it->first)
            {
                tag = it->second;
                break;
            }
       }

       if (tag == -1)
       {
            if (myRank == 0) std::cout << "No geometric group for boundary " << name << std::endl;
            continue;
       }

       type = read<ryml::csubstr>(settings,rymlKeys({"boundaryField",node.key(),"type"}));

       if (myRank == 0) std::cout << node.key() << "\t| " << tag << "\t| " << type << "; ";

       // set bdr markers
       bdr_markers[patchNumber] = 0;
       bdr_markers[patchNumber][tag-1] = 1;
       if (myRank == 0) std::cout << "set BC markers OK; ";

       // read additional type values
       if (type == "subsonicInletTotalPressure")
       {
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pTot"}));
            double inletT = read<double>(settings,rymlKeys({"boundaryField",node.key(),"TTot"}));

            A.AddBdrFaceIntegrator(new BoundaryIntegratorSubsonicInletTotalPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
       }
       else if (type == "subsonicInletFixedPressure")
       {
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pFix"}));
            double inletT = read<double>(settings,rymlKeys({"boundaryField",node.key(),"TTot"}));

            A.AddBdrFaceIntegrator(new BoundaryIntegratorSubsonicInletFixedPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
       }
       else if (type == "supersonicInlet")
       {
            double inletRho = read<double>(settings,rymlKeys({"boundaryField",node.key(),"rhoI"}));
            double inletU = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),0);
            double inletV = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),1);
            double inletW = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),2);
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pI"}));

            Vector inletVars(num_equation);

            inletVars(0) = inletRho;
            inletVars(1) = inletRho * inletU;
            inletVars(2) = inletRho * inletV;
            inletVars(3) = inletRho * inletW;
            inletVars(num_equation - 1) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;

            A.AddBdrFaceIntegrator(new BoundaryIntegratorSupersonicInlet(rsolver, dim, inletVars), bdr_markers[patchNumber] );
       }
       else if (type == "slip")
       {
            A.AddBdrFaceIntegrator(new BoundaryIntegratorSlip(rsolver, dim), bdr_markers[patchNumber] );
       }
       else if (type == "outlet")
       {
            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[patchNumber] );
       }
       else if (type == "totalPressureFixedTemperatureOutlet")
       {
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pTot"}));
            double inletT = read<double>(settings,rymlKeys({"boundaryField",node.key(),"TFix"}));

            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpenTotalPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
       }
       else if (type == "fixedPressureOutlet")
       {
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pFix"}));

            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpenFixedPressure(rsolver, dim, inletP), bdr_markers[patchNumber] );
       }
       else if (type == "charOutlet")
       {
            double inletRho = read<double>(settings,rymlKeys({"boundaryField",node.key(),"rhoI"}));
            double inletU = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),0);
            double inletV = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),1);
            double inletW = readVectorComponent<double>(settings,rymlKeys({"boundaryField",node.key(),"UI"}),2);
            double inletP = read<double>(settings,rymlKeys({"boundaryField",node.key(),"pI"}));

            Vector inletVars(num_equation);

            inletVars(0) = inletRho;
            inletVars(1) = inletRho * inletU;
            inletVars(2) = inletRho * inletV;
            inletVars(3) = inletRho * inletW;
            inletVars(num_equation - 1) = inletP / (specific_heat_ratio - 1.) + 0.5 * inletRho * (inletU * inletU + inletV * inletV)  ;

            A.AddBdrFaceIntegrator(new BoundaryIntegratorCharOutlet(rsolver, dim, inletVars), bdr_markers[patchNumber] );
       }
       else
       {
            if (myRank == 0)
            {
                std::cout << "Wrong boundary type " << type << std::endl;
                std::cout << "Available types: "
                        << "slip" << ", "
                        << "outlet" << ", "
                        << "subsonicInlet" << ", "
                        << "supersonicInlet" << ", "
                        << std::endl;
            } 
            exit(1);
       }
       if (myRank == 0) std::cout << "add boundary integrator OK" << std::endl;
       patchNumber++;
   }
}

void CaseManager::loadTimeSolver(ODESolver*& ode_solver, Limiter* l)
{
   int order = readOrDefault<int>(settings,rymlKeys({"time","order"}),0);

   if (myRank == 0) std::cout << "Runge --- Kutta order: " << order << std::endl;

   if (order == 2)
   {
       a = new double [3];
       a[0] = 0.0;
       a[1] = 1.0; 
       a[2] = 0.0;
       b = new double [2];
       b[0] = 0.5;
       b[1] = 0.5;
       c = new double [2];
       c[0] = 0.0;
       c[1] = 1.0;


       // alpha[2] = 1.0;

       // beta[0][0] = 0.3333333333333333;
            
       // beta[1][0] = 0.3333333333333333;
       // beta[1][1] = 0.3333333333333333;

       // beta[2][0] = 0.25;
       // beta[2][1] = 0.0;
       // beta[2][2] = 0.75;
   }
   else if (order == 3)
   {
       a = new double [6];
       a[0] = 0.0;
       // a[1] = 0.5; 
       // a[2] = 0.0;
       // a[3] = 0.0;
       // a[4] = 1.0; 
       // a[5] = 0.0;
       // b = new double [3];
       // b[0] = 0.1666666666666666;
       // b[1] = 0.6666666666666666;
       // b[2] = 0.1666666666666666;
       // c = new double [2];
       // c[0] = 0.0;
       // c[1] = 0.5;
       // c[2] = 1.0;
       a[0] = 0.0;
       a[1] = 1.0; 
       a[2] = 0.0;
       a[3] = 0.25;
       a[4] = 0.25; 
       a[5] = 0.0;
       b = new double [3];
       b[0] = 0.1666666666666666;
       b[1] = 0.1666666666666666;
       b[2] = 0.6666666666666666;
       c = new double [3];
       c[0] = 0.0;
       c[1] = 0.5;
       c[2] = 1.0;
   }
   else
   {
       if (myRank == 0) std::cout << "Wrong RK scheme order " << order << std::endl;
       exit(1);
   }

   ode_solver = new ExplicitRKLimitedSolver(order,a,b,c,*l);

   if (myRank == 0) std::cout << "RK time scheme order: " << order << std::endl;
}

void CaseManager::initializeRestartQueue()
{
   restart_current_cycle_name = "restart_000000";
   restart_deleted_cycle_name = "not_exists";

   if (restart) // change current name to current time cycle
   {
       restart_name_stream << std::setw(6) << std::setfill('0') << restart_data_c.GetCycle();
       restart_current_cycle_name = "restart_" + restart_name_stream.str();
       restart_name_stream.str(std::string());
   }
   
   if (myRank == 0)
   {
       // check old restart folders in the case directory and set them into the queue
       for (auto& p : std::filesystem::directory_iterator(caseDir))
            if (p.is_directory() 
                && 
                (p.path().filename().string().find("restart") != std::string::npos))
            {
                restart_deleted_cycle_name = p.path().filename().string();
                // check if some time cycles more than restarted
                if (restart_deleted_cycle_name.compare(restart_current_cycle_name) > 0)
                {
                     system(("rm -rf " + restart_deleted_cycle_name).c_str());
                     system(("rm -rf " + restart_deleted_cycle_name + ".mfem_root").c_str());
                }
                else
                {
                     restart_queue.push(restart_deleted_cycle_name);
                }
            }

       // std::cout << "--\n";
       // while (!restart_queue.empty())
       // {
       //    std::cout << restart_queue.top() << std::endl;
       //    restart_queue.pop();
       // }
   }
   // exit(0);

   if (restart_queue.size() == 0)
       restart_queue.push("restart_000000");
}

void CaseManager::loadTimeSettings(double& t_final, double& dt, double& cfl)
{
   t_final = read<double>(settings,rymlKeys({"time","tEnd"}));
   dt = read<double>(settings,rymlKeys({"time","dt"}));

   ryml::csubstr type = read<ryml::csubstr>(settings,rymlKeys({"time","type"}));

   if (type == "constant")
   {
       cfl = -1;
   }
   else
   {
       cfl = read<double>(settings,rymlKeys({"time","CoMax"}));
   }
}

double CaseManager::setStartTime()
{
   return restart ? restart_data_c.GetTime() : 0.0;
}

int CaseManager::setStartTimeCycle()
{
   return restart ? restart_data_c.GetCycle() : 0;
}

void CaseManager::cleanPreviousRestartFrames(int& ti)
{
   // // get current str name
   // restart_name_stream << std::setw(6) << std::setfill('0') << ti;
   // //std::cout << "restart_" + restart_name_stream.str() << std::endl;

   // restart_queue.push("restart_" + restart_name_stream.str());
   // if (restart_queue.size() > nSavedFrames)
   // {
   //    restart_deleted_cycle_name = restart_queue.top();
   //    system(("rm -rf " + restart_deleted_cycle_name).c_str());
   //    system(("rm -rf " + restart_deleted_cycle_name + ".mfem_root").c_str());
   //    restart_queue.pop();
   // }
   // restart_name_stream.str(std::string());
}

void CaseManager::getVisSteps(int& vis_steps)
{
   vis_steps = readOrDefault<int>(settings,rymlKeys({"postProcess","visSteps"}), 1);
   if (myRank == 0) std::cout << "vis_steps = " << vis_steps << std::endl;

   paraviewLevelOfDetails = readOrDefault<int>(settings,rymlKeys({"postProcess","levelOfDetails"}), 1);
   if (myRank == 0) std::cout << "level_of_details = " << paraviewLevelOfDetails << std::endl;
}