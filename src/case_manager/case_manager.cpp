#include <case_manager.hpp>

using namespace std;

CaseManager::CaseManager(std::string& caseFileName, VisItDataCollection& rdc)
:
    contents(NULL),
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

    // parse case file
    parse(caseFileName);

    // restart_queue
    // restart_name_stream
}

void CaseManager::parse(std::string& caseFileName)
{
    // read file
    ifstream reader;
    reader.open(caseFileName.c_str());
     
    if (!reader.is_open())
    {
        if (myRank == 0) cout << "Case file " << caseFileName << " is not found\n";
        exit(0);
    }
    else
    {
        if (myRank == 0)
        {
            cout << "===============================================" << endl;
            cout << "Reading case file " << caseFileName << "..." << endl;
            cout << "-----------------------------------------------" << endl;
        }
    };

    //parse result by ryml
    contents = new std::string((std::istreambuf_iterator<char>(reader)), 
        std::istreambuf_iterator<char>());

    reader.close();

    settings = new ryml::Tree(ryml::parse(c4::to_substr(*contents)));    

    // read restart settings
    c4::from_chars((*settings)["time"]["restart"].val(), &restart);
    c4::from_chars((*settings)["time"]["restartCycle"].val(), &restartCycle);
    c4::from_chars((*settings)["time"]["nSavedFrames"].val(), &nSavedFrames);
    if (myRank == 0) cout << "Restart of the case: " << restart << endl;

    // read physics
    std::string phType = ryml::preprocess_json<std::string>((*settings)["physics"]["type"].val());
    if (myRank == 0) cout << "Gas type: " << phType << endl;
    c4::from_chars((*settings)["physics"]["gamma"].val(), &specific_heat_ratio);
    if (myRank == 0) cout << "* Specific heat ratio: " << specific_heat_ratio << endl;
    if (phType == "covolume")
    {
        c4::from_chars((*settings)["physics"]["covolume"].val(), &covolume_constant);
        if (myRank == 0) cout << "* Covolume constant: " << covolume_constant << endl;
    }
    else
        covolume_constant = 0.0;

    
    if ((*settings)["physics"]["R"].val().size() != 0)
        c4::from_chars((*settings)["physics"]["R"].val(), &gas_constant);
    else
        gas_constant = 1.0;
    if (myRank == 0) cout << "* Gas constant: " << gas_constant << endl;
    
    // read dg settings
    c4::from_chars((*settings)["spatial"]["polyOrder"].val(), &spatialOrder);
    if (myRank == 0) cout << "Order of polynoms: " << spatialOrder << endl;

    //postprocess features
    
    if ((*settings)["postProcess"]["checkTotalEnergy"].val().size() != 0)
        c4::from_chars((*settings)["postProcess"]["checkTotalEnergy"].val(), &checkTotalEnergy);
    else
        checkTotalEnergy = 0;

    if ((*settings)["postProcess"]["writeIndicators"].val().size() != 0)
        c4::from_chars((*settings)["postProcess"]["writeIndicators"].val(), &writeIndicators);
    else
        writeIndicators = 0;


    //read initial conditions
    icType = ryml::preprocess_json<std::string>((*settings)["internalField"]["type"].val());
    if (myRank == 0) cout << "Type of problem: " << icType << endl;

    if (icType == "constant")
    {
        Vector sol(num_equation);
        c4::from_chars((*settings)["internalField"]["rho"].val(), &sol[0]);
        c4::from_chars((*settings)["internalField"]["U"][0].val(), &sol[1]);
        c4::from_chars((*settings)["internalField"]["U"][1].val(), &sol[2]);
        c4::from_chars((*settings)["internalField"]["U"][2].val(), &sol[3]);
        c4::from_chars((*settings)["internalField"]["p"].val(), &sol[num_equation-1]);

        ICInterface = new ICConstant(sol);
    }
    else if (icType == "planeBreakup")
    {
        Vector sol1(num_equation);
        Vector sol2(num_equation);
        

        c4::from_chars((*settings)["internalField"]["origin"][0].val(), &origin[0]);
        c4::from_chars((*settings)["internalField"]["origin"][1].val(), &origin[1]);
        c4::from_chars((*settings)["internalField"]["origin"][2].val(), &origin[2]);
        c4::from_chars((*settings)["internalField"]["normal"][0].val(), &normal[0]);
        c4::from_chars((*settings)["internalField"]["normal"][1].val(), &normal[1]);
        c4::from_chars((*settings)["internalField"]["normal"][2].val(), &normal[2]);

        c4::from_chars((*settings)["internalField"]["left"]["rho"].val(), &sol1[0]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][0].val(), &sol1[1]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][1].val(), &sol1[2]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][2].val(), &sol1[3]);
        c4::from_chars((*settings)["internalField"]["left"]["p"].val(), &sol1[num_equation-1]);

        c4::from_chars((*settings)["internalField"]["right"]["rho"].val(), &sol2[0]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][0].val(), &sol2[1]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][1].val(), &sol2[2]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][2].val(), &sol2[3]);
        c4::from_chars((*settings)["internalField"]["right"]["p"].val(), &sol2[num_equation-1]);

        ICInterface = new ICPlaneBreakup(sol1, sol2, origin, normal);
    }
    else if (icType == "sphericalBreakup")
    {
        Vector sol1(num_equation);
        Vector sol2(num_equation);
        double radius = 0.0;

        c4::from_chars((*settings)["internalField"]["origin"][0].val(), &origin[0]);
        c4::from_chars((*settings)["internalField"]["origin"][1].val(), &origin[1]);
        c4::from_chars((*settings)["internalField"]["origin"][2].val(), &origin[2]);
        c4::from_chars((*settings)["internalField"]["radius"].val(), &radius);

        c4::from_chars((*settings)["internalField"]["left"]["rho"].val(), &sol1[0]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][0].val(), &sol1[1]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][1].val(), &sol1[2]);
        c4::from_chars((*settings)["internalField"]["left"]["U"][2].val(), &sol1[3]);
        c4::from_chars((*settings)["internalField"]["left"]["p"].val(), &sol1[num_equation-1]);

        c4::from_chars((*settings)["internalField"]["right"]["rho"].val(), &sol2[0]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][0].val(), &sol2[1]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][1].val(), &sol2[2]);
        c4::from_chars((*settings)["internalField"]["right"]["U"][2].val(), &sol2[3]);
        c4::from_chars((*settings)["internalField"]["right"]["p"].val(), &sol2[num_equation-1]);

        ICInterface = new ICSphericalBreakup(sol1, sol2, origin, radius);
    }
    else if (icType == "densityPulse")
    {
        Vector sol(num_equation);
        double epsilon = 0.0;

        c4::from_chars((*settings)["internalField"]["origin"][0].val(), &origin[0]);
        c4::from_chars((*settings)["internalField"]["origin"][1].val(), &origin[1]);
        c4::from_chars((*settings)["internalField"]["origin"][2].val(), &origin[2]);
        c4::from_chars((*settings)["internalField"]["epsilon"].val(), &epsilon);

        c4::from_chars((*settings)["internalField"]["U"][0].val(), &sol[1]);
        c4::from_chars((*settings)["internalField"]["U"][1].val(), &sol[2]);
        c4::from_chars((*settings)["internalField"]["U"][2].val(), &sol[3]);
        c4::from_chars((*settings)["internalField"]["p"].val(), &sol[num_equation-1]);

        ICInterface = new ICDensityPulse(sol, origin, epsilon);
    }
    else
    {
        if (myRank == 0)
        {
            cout << "Wrong type if initial condition" << endl
             << "Available types: "
             << "constant, planeBreakup, sphericalBreakup" << endl;
        }
    }

    if (myRank == 0) cout << "Initial conditions OK" << endl;
}

void CaseManager::loadMesh(ParMesh*& pmesh)
{
    const char* meshFile = ryml::preprocess_json<std::string>((*settings)["mesh"]["file"].val()).c_str();
    
    // read mesh settings
    // strcpy(meshFile, mf);
    c4::from_chars((*settings)["mesh"]["serialRefLevels"].val(), &serRefLevels);
    c4::from_chars((*settings)["mesh"]["parallelRefLevels"].val(), &parRefLevels);
    c4::from_chars((*settings)["mesh"]["adaptive"].val(), &adaptive_mesh);

    if (myRank == 0)
    {
        cout << "Reading mesh from " << ryml::preprocess_json<std::string>((*settings)["mesh"]["file"].val()).c_str() << "..." << endl;
        cout << "* Serial refinement levels = " << serRefLevels << endl;
        cout << "* Parallel refinement levels = " << parRefLevels << endl;
        cout << "* Adaptive: " << (adaptive_mesh ? "true" : "false") << endl;
    }
    // read adaptive mesh settings


    if (restart) // load par mesh in case of restart
    {
        // 2.1. Read the given data collection.
        restart_data_c.Load(restartCycle);
        if (restart_data_c.Error())
        {
            if (myRank == 0)
            {
                cout << "Error loading data collection" << endl;
            }
            return;
        }
        // 2.2. Get par mesh pieces from the loaded collection
        pmesh = dynamic_cast<ParMesh*>(restart_data_c.GetMesh());
        if (!pmesh)
        {
            if (myRank == 0)
            {
                cout << "The given data collection does not have a parallel mesh."
                      << endl;
            }
            return;
        }
    }
    else //load mesh from file and make a partition
    {
        // 2.1. Read the serial mesh on all processors, refine it in serial, then
        //     partition it across all processors and refine it in parallel.

        std::ifstream meshStream(ryml::preprocess_json<std::string>((*settings)["mesh"]["file"].val()).c_str());

        Mesh *mesh = new Mesh(meshStream, 1, 1);
        mesh->EnsureNCMesh(true);

        for (int l = 0; l < serRefLevels; l++)
        {
            mesh->UniformRefinement();
        }

        //if (myRank == 0) { cout << "Number of cells: " << mesh->GetNE() << endl; }

        // 2.2. Define a parallel mesh by a partitioning of the serial mesh. Refine
        //     this mesh further in parallel to increase the resolution. Once the
        //     parallel mesh is defined, the serial mesh can be deleted.
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;

        for (int l = 0; l < parRefLevels; l++)
        {
            pmesh->UniformRefinement();
        }
        restart_data_c.SetFormat(1);
        restart_data_c.SetMesh(pmesh); 
        if (myRank == 0) cout << "Number of cells: " << pmesh->GetNE() << endl; 
    } 

    if (myRank == 0) cout << pmesh->Dimension() << endl;

    readPhysicalNames(pmesh);
    minimizeAttributes(pmesh);
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

    std::ifstream input (ryml::preprocess_json<std::string>((*settings)["mesh"]["file"].val()).c_str());

    while (input >> buff)
    {
        if (buff == "$PhysicalNames") // reading mesh vertices
        {
            input >> number_of_physical_groups;
            getline(input, buff);

            for (int iGroup = 0; iGroup < number_of_physical_groups; ++iGroup)
            {
                input >> dim_r >> tag_r >> phys_group_name_r;

                if ( phys_group_name_r.front() == '"' ) {
                    phys_group_name_r.erase( 0, 1 ); // erase the first character
                    phys_group_name_r.erase( phys_group_name_r.size() - 1 ); // erase the last character
                }

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
                    cout << "Wrong matching of dimension for physical group"
                         << "(" << phys_group_name_r << ":" << dim_r << ") "
                         << "and space dimension " 
                         << mesh->SpaceDimension()
                         << endl;
                    exit(1);
                }
            }

            break;
        }
        else
        if (buff == "$EndElements")
        {
            cout << "No physical boundaries in this mesh";
            break;
        }
    }

    input.close();
}

void CaseManager::minimizeAttributes(ParMesh*& mesh)
{
    // start from internal cells; default attribute should be 1
    mesh->bdr_attributes.Print(cout << "mesh.bdr_attributes = ");
    mesh->attributes.Print(cout << "mesh.attributes = ");

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
            //     hasUnitBdrAttribute = 1;
            //     break;
            // }
            if (curAttr == mesh->bdr_attributes[iAttr])
            {
                mesh->SetBdrAttribute(iCell, iAttr + 1 + mesh->attributes.Size());
                break;
            }
        }
    }

    // cout << "hasUnitBdrAttribute = " << hasUnitBdrAttribute << endl;

    // edit boundary attributes inside map for names
    std::map<std::string,int>& map_bdr_names_tag = (mesh->SpaceDimension() == 2) ? map_phys_names_tag_1D : map_phys_names_tag_2D;

    for(std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
    {
        for (int iAttr = 0; iAttr < mesh->bdr_attributes.Size(); ++iAttr)
        {
            if (it->second == mesh->bdr_attributes[iAttr])
            {
                it->second = iAttr + 1 + mesh->attributes.Size();
                break;
            }
        }
    }

    // edit mesh attributes
    // for (int i = 0; i < mesh->bdr_attributes.Size() - hasUnitBdrAttribute; ++i)
    //     mesh->bdr_attributes[i + hasUnitBdrAttribute] = i + 1 + mesh->attributes.Size();
    
    for (int i = 0; i < mesh->bdr_attributes.Size(); ++i)
        mesh->bdr_attributes[i] = i + 1 + mesh->attributes.Size();

    // if (hasUnitBdrAttribute)
    //     mesh->bdr_attributes[0] = 1;

}

void CaseManager::loadAdaptiveMeshSettings(ThresholdRefiner& refiner,ThresholdDerefiner& derefiner)
{
    c4::from_chars((*settings)["mesh"]["totalErrorFraction"].val(), &total_error_fraction);
    c4::from_chars((*settings)["mesh"]["maxElemError"].val(), &max_elem_error);
    c4::from_chars((*settings)["mesh"]["hysteresis"].val(), &hysteresis);
    c4::from_chars((*settings)["mesh"]["preferConformingRefinement"].val(), &prefer_conforming_refinement);
    c4::from_chars((*settings)["mesh"]["nonConformingLimit"].val(), &nc_limit);

    refiner.SetTotalErrorFraction(total_error_fraction); // use purely local threshold

    refiner.SetLocalErrorGoal(max_elem_error);
    if (prefer_conforming_refinement)
        refiner.PreferConformingRefinement();
    refiner.SetNCLimit(nc_limit);

    derefiner.SetThreshold(hysteresis * max_elem_error);
    derefiner.SetNCLimit(nc_limit);
}

void CaseManager::checkNumEqn(int dim)
{
    if (dim == 2)
        num_equation = 4;
    else if (dim == 3)
        num_equation = 5;
    else
    {
        if (myRank == 0)
        {
            cout << "Wrong dimension " << dim << endl;
            cout << "Please check your mesh file" << endl;
        }
        exit(1);
    }

    if (myRank == 0) cout << "Number of equations = " << num_equation << endl;
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
                cout << "The given data collection has no 'restart_conservative_vars' field."
                      << endl;
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
        Vector solConst(num_equation);
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
    rsolverType = ryml::preprocess_json<std::string>((*settings)["spatial"]["riemannSolver"].val());

    if (myRank == 0) cout << "Riemann solver type: " << rsolverType << endl;

    if (rsolverType == "Rusanov")
        rsolver = new RiemannSolverRusanov();
    else if (rsolverType == "LLF")
        rsolver = new RiemannSolverLLF();
    else if (rsolverType == "HLL")
        rsolver = new RiemannSolverHLL();
    else if (rsolverType == "HLLC")
        rsolver = new RiemannSolverHLLC();
    else 
    {
        if (myRank == 0) cout << "Unknown Riemann solver type: " << rsolverType << '\n';
        exit(1);
    }
}

void CaseManager::loadLimiter(Averager& avgr, Indicator*& ind, Limiter*& l, const Array<int>& offsets, const int dim, BlockVector& indicatorData, ParFiniteElementSpace& vfes)
{
    std::string limiterType = ryml::preprocess_json<std::string>((*settings)["spatial"]["limiter"]["type"].val());

    std::string indicatorType = "notype";

    if ((*settings)["spatial"]["limiter"]["linearize"].val().size() != 0)
        c4::from_chars((*settings)["spatial"]["limiter"]["linearize"].val(), &linearize);
    else
        linearize = 0;

    if ((*settings)["spatial"]["limiter"]["haveLastHope"].val().size() != 0)
        c4::from_chars((*settings)["spatial"]["limiter"]["haveLastHope"].val(), &haveLastHope);
    else
        haveLastHope = 1;

    cout << "Additional linearization: " << linearize << endl;
    cout << "Cut slopes in case of nonphysical values in vertices: " << haveLastHope << endl;
    // if (limiterType == "BJ") 
    // {
    //     indicatorType = "BJ";
    //     limiterType = "Multi";
    // }
    // else
    // {
        indicatorType = ryml::preprocess_json<std::string>((*settings)["spatial"]["limiter"]["indicator"].val());
    // }

    if (indicatorType == "Nowhere")
    {
       ind = new IndicatorNowhere(avgr, &vfes, offsets, dim, indicatorData);
    }
    else if (indicatorType == "Everywhere")
    {
        ind = new IndicatorEverywhere(avgr, &vfes, offsets, dim, indicatorData); 
    }
    else if (indicatorType == "BJ")
    {
        ind = new IndicatorBJ(avgr, &vfes, offsets, dim, indicatorData); 
    }
    else if (indicatorType == "Shu")
    {
        ind = new IndicatorShu(avgr, &vfes, offsets, dim, indicatorData);    
    }
    else
    {
        if (myRank == 0) cout << "Unknown indicator type: " << indicatorType << '\n';
        exit(1);
    }  

    // Define limiter object
    
    if (limiterType == "FinDiff")
    {
        l = new LimiterFinDiff(*ind, avgr, &vfes, offsets, linearize, haveLastHope, dim);
    }
    else if (limiterType == "Multi")
    {
        l = new LimiterMultiplier(*ind, avgr, &vfes, offsets, linearize, haveLastHope, dim);
    }
    else
    {
        if (myRank == 0) cout << "Unknown limiter type: " << limiterType << '\n';
        exit(1);
    }

    if (myRank == 0) cout << "Indicator type: " << indicatorType << "; Limiter type: " << limiterType << endl;
}

void CaseManager::addBoundaryIntegrators(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim)
{
    int tag = 0;
    std::string name;
    std::string type;
     
    // for fully periodic meshes
    if (mesh.bdr_attributes.Size() == 0) 
        return;

    std::map<std::string,int>& map_phys_names_tag = (mesh.SpaceDimension() == 2) ? map_phys_names_tag_2D : map_phys_names_tag_3D;
    std::map<std::string,int>& map_bdr_names_tag  = (mesh.SpaceDimension() == 2) ? map_phys_names_tag_1D : map_phys_names_tag_2D;

    if (myRank == 0) 
    {
        mesh.bdr_attributes.Print(cout << "mesh.bdr_attributes = ");
        mesh.attributes.Print(cout << "mesh.attributes = ");

        cout << "Boundary mapping:\n";
        for ( std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
        {
            cout << it->first << ' ' << it->second << endl;
        }
        cout << "Physical group mapping:\n";
        for ( std::map<std::string,int>::iterator it = map_phys_names_tag.begin(); it != map_phys_names_tag.end(); it++)
        {
            cout << it->first << ' ' << it->second << endl;
        } 
    }


    bdr_markers.resize(mesh.bdr_attributes.Max());
    for (int i = 0; i < mesh.bdr_attributes.Max(); ++i)
        bdr_markers[i].SetSize(mesh.bdr_attributes.Max());

    // if (myRank == 0) mesh.bdr_attributes.Print(cout << "Boundary attributes in mesh:\n");   

    

    // if ( ((*settings)["boundaryField"].num_children() != mesh.bdr_attributes.Size()) )
    // {
    //     if (myRank == 0)
    //     {
    //         cout << "=====\nPossible mismatch of boundary conditions and geometric boundary groups." << endl;
    //         cout << "(*settings)[boundaryField].num_children() = " << (*settings)["boundaryField"].num_children() << endl;
    //         cout << "mesh.bdr_attributes.Size() = " << mesh.bdr_attributes.Size() << endl;
    //         cout << "Boundary names:" << endl;
    //         for (ryml::NodeRef node : (*settings)["boundaryField"].children())
    //         {
    //             cout << " * " << node.key() << endl;
    //         }
    //         cout << "Geometric boundaries:\n";
    //         for ( std::map<std::string,int>::iterator it = map_bdr_names_tag.begin(); it != map_bdr_names_tag.end(); it++)
    //         {
    //             cout << " * " << it->first << endl;
    //         }
    //     }
        
    //     exit(1);
    // }

    if (myRank == 0) cout << "Boundary conditions:" << endl;

    // read patch names and tags from mesh file $PhysicalNames
    int patchNumber = 0;
    for (ryml::NodeRef node : (*settings)["boundaryField"].children())
    {
        // find tag for the name
        tag = -1;
        name = ryml::preprocess_json<std::string>(node.key());
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
            if (myRank == 0) cout << "No geometric group for boundary " << name << endl;
            continue;
        }
        // read names, tags and sth like this
        c4::from_chars((*settings)["boundaryField"][node.key()]["type"].val(), &type);

        if (myRank == 0) cout << node.key() << "\t| " << tag << "\t| " << type << "; ";

        // set bdr markers
        bdr_markers[patchNumber] = 0;
        bdr_markers[patchNumber][tag-1] = 1;
        if (myRank == 0) cout << "set BC markers OK; ";

        // read additional type values
        if (type == "subsonicInletTotalPressure")
        {
            double inletT;
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["pTot"].val(), &inletP);
            c4::from_chars((*settings)["boundaryField"][node.key()]["TTot"].val(), &inletT);

            A.AddBdrFaceIntegrator(new BoundaryIntegratorSubsonicInletTotalPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
        }
        else if (type == "subsonicInletFixedPressure")
        {
            double inletT;
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["pFix"].val(), &inletP);
            c4::from_chars((*settings)["boundaryField"][node.key()]["TTot"].val(), &inletT);

            A.AddBdrFaceIntegrator(new BoundaryIntegratorSubsonicInletFixedPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
        }
        else if (type == "supersonicInlet")
        {
            double inletRho;
            double inletU;
            double inletV;
            double inletW;
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["rhoI"].val(), &inletRho);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][0].val(), &inletU);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][1].val(), &inletV);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][2].val(), &inletW);
            c4::from_chars((*settings)["boundaryField"][node.key()]["pI"].val(), &inletP);

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
            double inletT;
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["pTot"].val(), &inletP);
            c4::from_chars((*settings)["boundaryField"][node.key()]["TFix"].val(), &inletT);

            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpenTotalPressure(rsolver, dim, inletP, inletT), bdr_markers[patchNumber] );
        }
        else if (type == "fixedPressureOutlet")
        {
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["pFix"].val(), &inletP);

            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpenFixedPressure(rsolver, dim, inletP), bdr_markers[patchNumber] );
        }
        else if (type == "charOutlet")
        {
            double inletRho;
            double inletU;
            double inletV;
            double inletW;
            double inletP;

            c4::from_chars((*settings)["boundaryField"][node.key()]["rhoI"].val(), &inletRho);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][0].val(), &inletU);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][1].val(), &inletV);
            c4::from_chars((*settings)["boundaryField"][node.key()]["UI"][2].val(), &inletW);
            c4::from_chars((*settings)["boundaryField"][node.key()]["pI"].val(), &inletP);

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
                std::cout << "Wrong boundary type " << type << endl;
                std::cout << "Available types: "
                      << "slip" << ", "
                      << "outlet" << ", "
                      << "subsonicInlet" << ", "
                      << "supersonicInlet" << ", "
                      << endl;
            } 
            exit(1);
        }
        if (myRank == 0) cout << "add boundary integrator OK" << endl;
        patchNumber++;
    }
}

void CaseManager::loadTimeSolver(ODESolver*& ode_solver, Limiter* l)
{
    int order = 0;
    c4::from_chars((*settings)["time"]["order"].val(), &order);

    if (myRank == 0) cout << "Runge --- Kutta order: " << order << endl;

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
        c = new double [2];
        c[0] = 0.0;
        c[1] = 0.5;
        c[2] = 1.0;
    }
    else
    {
        if (myRank == 0) cout << "Wrong RK scheme order " << order << endl;
        exit(1);
    }

    ode_solver = new ExplicitRKLimitedSolver(order,a,b,c,*l);

    if (myRank == 0) cout << "RK time scheme order: " << order << endl;
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

        // cout << "--\n";
        // while (!restart_queue.empty())
        // {
        //     cout << restart_queue.top() << endl;
        //     restart_queue.pop();
        // }
    }
    // exit(0);

    if (restart_queue.size() == 0)
        restart_queue.push("restart_000000");
}

void CaseManager::loadTimeSettings(double& t_final, double& dt, double& cfl)
{
    c4::from_chars((*settings)["time"]["tEnd"].val(), &t_final);
    c4::from_chars((*settings)["time"]["dt"].val(), &dt);

    std::string type = ryml::preprocess_json<std::string>((*settings)["time"]["type"].val());

    if (type == "constant")
    {
        cfl = -1;
    }
    else
    {
        c4::from_chars((*settings)["time"]["CoMax"].val(), &cfl);
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
    // get current str name
    restart_name_stream << std::setw(6) << std::setfill('0') << ti;
    //cout << "restart_" + restart_name_stream.str() << endl;

    restart_queue.push("restart_" + restart_name_stream.str());
    if (restart_queue.size() > nSavedFrames)
    {
        restart_deleted_cycle_name = restart_queue.top();
        system(("rm -rf " + restart_deleted_cycle_name).c_str());
        system(("rm -rf " + restart_deleted_cycle_name + ".mfem_root").c_str());
        restart_queue.pop();
    }
    restart_name_stream.str(std::string());
}

void CaseManager::getVisSteps(int& vis_steps)
{
    c4::from_chars((*settings)["postProcess"]["visSteps"].val(), &vis_steps);

    if (myRank == 0) cout << "vis_steps = " << vis_steps << endl;
}