#include <case_manager.hpp>

using namespace std;

CaseManager::CaseManager(std::string& caseFileName, VisItDataCollection& rdc)
:
    restart_data_c (rdc),
    origin(num_equation-2),
    normal(num_equation-2)
{
    caseDir = std::filesystem::current_path().string();

    // set default values for parameters
    restart = false;
    restartCycle = 0;
    nSavedFrames = 0;
    icType = "constant";
    rsolverType = "Rusanov";

    caseFileName = caseDir + "/" + caseFileName;

    total_error_fraction = 0.0;
    max_elem_error = 0.0;
    hysteresis = 0.0;
    nc_limit  = 0;
    prefer_conforming_refinement = false;

    // parse case file
    parse(caseFileName);
}

void CaseManager::parse(std::string& caseFileName)
{
    // read file
    ifstream reader;
    reader.open(caseFileName.c_str());
     
    if (!reader.is_open())
    {
        cout << "Case file " << caseFileName << " is not found\n";
        exit(0);
    }
    else
    {
        cout << "===============================================" << endl;
        cout << "Reading case file " << caseFileName << "..." << endl;
        cout << "-----------------------------------------------" << endl;
    };

    //parse result by ryml
    contents = new std::string((std::istreambuf_iterator<char>(reader)), 
        std::istreambuf_iterator<char>());

    reader.close();

    settings = new ryml::Tree(ryml::parse(c4::to_substr(*contents)));    

    // read restart settings
    c4::from_chars((*settings)["time"]["restart"].val(), &restart);
    c4::from_chars((*settings)["time"]["restartCycle"].val(), &restartCycle);
    cout << "Restart case: " << restart << endl;

    // read physics
    std::string phType = ryml::preprocess_json<std::string>((*settings)["physics"]["type"].val());
    cout << "Gas type: " << phType << endl;
    c4::from_chars((*settings)["physics"]["gamma"].val(), &specific_heat_ratio);
    cout << "* Specific heat ratio: " << specific_heat_ratio << endl;
    if (phType == "covolume")
    {
        c4::from_chars((*settings)["physics"]["covolume"].val(), &covolume_constant);
        cout << "* Covolume constant: " << covolume_constant << endl;
    }
    else
        covolume_constant = 0.0;
    
    // read dg settings
    c4::from_chars((*settings)["spatial"]["polyOrder"].val(), &spatialOrder);
    cout << "Order of polynoms: " << spatialOrder << endl;


    //read initial conditions
    icType = ryml::preprocess_json<std::string>((*settings)["internalField"]["type"].val());
    cout << "Type of problem: " << icType << endl;

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
    else
    {
        cout << "Wrong type if initial condition" << endl
             << "Available types: "
             << "constant, planeBreakup, sphericalBreakup" << endl;
    }

    cout << "Initial conditions OK" << endl;
}

void CaseManager::loadMesh(ParMesh*& pmesh)
{
    const char* meshFile = ryml::preprocess_json<std::string>((*settings)["mesh"]["file"].val()).c_str();
    
    // read mesh settings
    // strcpy(meshFile, mf);
    c4::from_chars((*settings)["mesh"]["serialRefLevels"].val(), &serRefLevels);
    c4::from_chars((*settings)["mesh"]["parallelRefLevels"].val(), &parRefLevels);
    c4::from_chars((*settings)["mesh"]["adaptive"].val(), &adaptive_mesh);

    cout << "Reading mesh from " << meshFile << "..." << endl;
    cout << "* Serial refinement levels = " << serRefLevels << endl;
    cout << "* Parallel refinement levels = " << parRefLevels << endl;
    cout << "* Adaptive: " << (adaptive_mesh ? "true" : "false") << endl;

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

        Mesh *mesh = new Mesh(meshFile, 1, 1);
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
        cout << "Number of cells: " << pmesh->GetNE() << endl; 
    } 

    cout << pmesh->Dimension() << endl;
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
        cout << "Wrong dimension " << dim << endl;
        cout << "Please check your mesh file" << endl;
        exit(1);
    }

    cout << "Number of equations = " << num_equation << endl;
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
        Vector solConst(num_equation);
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
    rsolverType = ryml::preprocess_json<std::string>((*settings)["spatial"]["riemannSolver"].val());

    cout << "Riemann solver type: " << rsolverType << endl;

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
        cout << "Unknown Riemann solver type: " << rsolverType << '\n';
        exit(1);
    }
}

void CaseManager::loadLimiter(Averager& avgr, Indicator*& ind, Limiter*& l, const Array<int>& offsets, const int dim, BlockVector& indicatorData, ParFiniteElementSpace& vfes)
{
    std::string limiterType = ryml::preprocess_json<std::string>((*settings)["spatial"]["limiter"]["type"].val());

    std::string indicatorType = "notype";
    if (limiterType == "BJ") 
    {
        indicatorType = "BJ";
        limiterType = "Multi";
    }
    else
    {
        indicatorType = ryml::preprocess_json<std::string>((*settings)["spatial"]["limiter"]["indicator"].val());
    }

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
        cout << "Unknown indicator type: " << indicatorType << '\n';
        exit(1);
    }  

    // Define limiter object
    
    if (limiterType == "FinDiff")
    {
        l = new LimiterFinDiff(*ind, avgr, &vfes, offsets, dim);
    }
    else if (limiterType == "Multi")
    {
        l = new LimiterMultiplier(*ind, avgr, &vfes, offsets, dim);
    }
    else
    {
        cout << "Unknown limiter type: " << limiterType << '\n';
        exit(1);
    }

    cout << "Indicator type: " << indicatorType << "; Limiter type: " << limiterType << endl;
}

void CaseManager::addBoundaryIntegrators(ParNonlinearForm& A, ParMesh& mesh, RiemannSolver& rsolver, int dim)
{
    int tag = 0;
    std::string type;

    bdr_markers.resize(mesh.bdr_attributes.Max());
    for (int i = 0; i < mesh.bdr_attributes.Max(); ++i)
        bdr_markers[i].SetSize(mesh.bdr_attributes.Max());

    cout << "Boundary conditions:" << endl;

    // read patch names and tags from mesh file $PhysicalNames
    int patchNumber = 0;
    for (ryml::NodeRef node : (*settings)["boundaryField"].children())
    {
        // read names, tags and sth like this
        c4::from_chars((*settings)["boundaryField"][node.key()]["tag"].val(), &tag);
        c4::from_chars((*settings)["boundaryField"][node.key()]["type"].val(), &type);

        cout << node.key() << "\t|" << type << endl;
        // set bdr markers
        

        bdr_markers[patchNumber] = 0;
        

        bdr_markers[patchNumber][tag] = 1;
        cout << "Set BC markers OK" << endl;


        // read additional type values
        if (type == "inlet")
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

            A.AddBdrFaceIntegrator(new BoundaryIntegratorConstant(rsolver, dim, inletVars), bdr_markers[patchNumber] );
        }
        else if (type == "slip")
        {
            A.AddBdrFaceIntegrator(new BoundaryIntegratorWall(rsolver, dim), bdr_markers[patchNumber] );
        }
        else if (type == "outlet")
        {
            A.AddBdrFaceIntegrator(new BoundaryIntegratorOpen(rsolver, dim), bdr_markers[patchNumber] );
        }
        else
        {
            std::cout << "Wrong boundary type " << type << endl;
            exit(1);
        }
        cout << "Add BI OK" << endl;
        patchNumber++;
    }
}

void CaseManager::loadTimeSolver(ODESolver*& ode_solver, Limiter* l)
{
    int order = 0;
    c4::from_chars((*settings)["time"]["order"].val(), &order);

    cout << "Runge --- Kutta order: " << order << endl;

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
        cout << "Wrong RK scheme order " << order << endl;
        exit(1);
    }

    ode_solver = new ExplicitRKLimitedSolver(order,a,b,c,*l);

    cout << "RK time scheme order: " << order << endl;
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

    cout << "vis_steps = " << vis_steps << endl;
}