Info<< "\nLooking for FLmap flag in elmerToFoamDict\n" << endl;

IOdictionary elmerToFoamDict
(
    IOobject
    (
        "elmerToFoamDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// These variables have to be declared outside the if(FLmap_) for Ueqn.H
const bool QJmap_ = elmerToFoamDict.getOrDefault<label>("QJmap", 0);
const bool QJcoef_ = elmerToFoamDict.getOrDefault<label>("QJcoef", 0);
const bool FLmap_ = elmerToFoamDict.getOrDefault<label>("FLmap", 0);
const scalar FLcoef_ = elmerToFoamDict.getOrDefault<label>("FLcoef", 1);

// Settings for rotation of 2D Elmer data to 3D OpenFOAM data
// should QJ and FL be rotated
const bool rot2D_ = elmerToFoamDict.getOrDefault<label>("rot2D", 0);
// axis of rotation in OpenFOAM
const vector axis_ = elmerToFoamDict.getOrDefault<vector>("axis",vector(0,1,0));
// direction of rotation in OpenFOAM
const vector rotDir_ = elmerToFoamDict.getOrDefault<vector>("rotDir",vector(0,0,1));

dimensionedScalar  FLrho0_
(
    elmerToFoamDict.getOrDefault<dimensionedScalar>
    (
        "FLrho0",
        1
    )
);
dimensionedScalar  Qrhocp_
(
    elmerToFoamDict.getOrDefault<dimensionedScalar>
    (
        "Qrhocp",
        1
    )
);

Info<< "  QJmap = " << QJmap_ << endl;
Info<< "  QJcoef = " << QJcoef_ << endl;
Info<< "  FLrho0 = " << FLrho0_ << endl;
Info<< "  FLmap = " << FLmap_ << endl;
Info<< "  FLcoef = " << FLcoef_ << endl;
Info<< "  FLrho0 = " << FLrho0_ << endl;

Info<< "Creating field QJ\n" << endl;

volScalarField QJ
(
    IOobject
    (
        "QJ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("dummyvector1", dimensionSet(1,-1,-3,0,0,0,0), 0 )
);

Info<< "Creating field FL\n" << endl;

volVectorField FL
(
    IOobject
    (
        "FL",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector ("dummyvector1", dimensionSet(1,-2,-2,0,0,0,0), vector(0,0,0) )
);

if (FLmap_ || QJmap_)
{
    const label forceColumn_(elmerToFoamDict.get<label>("forceColumn"));
    const label heatColumn_(elmerToFoamDict.get<label>("heatColumn"));
    const label coordColumn_(elmerToFoamDict.get<label>("coordColumn"));
    const word FLfile_ = elmerToFoamDict.get<word>("FLfile");
    const word method_ = elmerToFoamDict.getOrDefault<word>("mapMethod", "interp1");

    // get labels for field rotation from axis_ and rotDir_
    label rrLabel;
    label rotLabel;
    label axisLabel;
    if (rot2D_)
    {
        rrLabel = (vector(1,1,1)-axis_-rotDir_) & vector(0,1,2);
        axisLabel = axis_ & vector(0,1,2);    
        rotLabel = rotDir_ & vector(0,1,2);
        Info << rrLabel << ' '<< axisLabel << ' '<< rotLabel << ' ' << endl;
    }

    fileName posName(FLfile_);
    IFstream inFile(posName);
    
    if(QJmap_)
    {
        DynamicList<vector> coordinates;
        DynamicList<scalar> heat;
        Info << "  Reading file "  << FLfile_ << endl;
        Info << "  coordinates at columns (x,y,z) = ["
            << coordColumn_ << ", " << coordColumn_+1 << ", " << coordColumn_+2 
            << "]" << endl;
        Info << "  Joule heat at column (QJ) = ["
            << heatColumn_ << "]" 
            << endl;
            
        while (inFile.good())
        {
            string line;
            inFile.getLine(line);

            if (line.size()>1)
            {
                // Set up up a stream from this line
                std::stringstream lineStr(line);

                DynamicList<scalar> data;
                string item;
                // Info << "line = " << line << endl;

                std::string token;
                while (std::getline(lineStr, token, ',')) 
                {            
                    double value = std::stod(token);  // Convert the token to a double
                    data.append(value);      // Store the value in the vector
                }

                // read Elmer coordinates
                //TBD: add permutation
                vector coord(Zero);
                forAll (coord, xi)
                {
                    coord[xi] = data[coordColumn_+xi];
                }   
                
                coordinates.append(coord);
                heat.append(data[heatColumn_]);
            }
        }
        inFile.rewind();

        Info<< "  coordinates.size() = " << coordinates.size() << endl;
        Info<< "  heat.size() = " << heat.size() << endl;
        Info<< "  Interpolation method: " << method_ << endl;
        
        // Mapping
        Info<< "\nMapping data\n" << endl;

        scalar maxdistall = -1000;

        forAll(mesh.cells(), cellsi)
        {
            vector p = mesh.C()[cellsi];
            vector por = p;
            // calculte r and project mesh points to 2D
            scalar rr = 0;
            // calculate radius rr from axis_
            if(rot2D_)
            {
                forAll (p, xi)
                {
                    rr = rr + (1-axis_[xi])*p[xi]*p[xi];
                }
                rr = Foam::sqrt(rr);
                // projected p
                p[rrLabel] = rr;
                p[axisLabel] = p[axisLabel];
                p[rotLabel] = 0;
            }

            // find nearest neighbour
            scalar mindist = 1000;
            label minpoint = 0;
            for (label i=0; i<coordinates.size(); i++) 
            {
                scalar dist = mag(p-coordinates[i]);
                if (dist < mindist) 
                { 
                    mindist = dist;  
                    minpoint = i; 
                }
            }
            // Info << "por = " << por << "; p = " << p << "; minpoint = " << coordinates[minpoint] << endl; 
            if (mindist>maxdistall) 
            { 
                maxdistall = mindist; 
            }     

            scalar interpHeat(Zero);

            if (method_ =="nearest")
            {
                interpHeat = heat[minpoint];
            }
            else if (method_ == "interp1")
            {
                // find 3 nearest points to p
                label minpoint1=0, minpoint2=0, minpoint3=0;
                // point1
                scalar mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) 
                    { 
                        mindist = dist;  
                        minpoint1 = i; 
                    }
                }
                if (mindist>maxdistall) { maxdistall = mindist; }
                // point2
                mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist && i!=minpoint1) 
                        { 
                            mindist = dist;  
                            minpoint2 = i; 
                        }
                    }
                if (mindist>maxdistall) { maxdistall = mindist; }

                // point3
                mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist && i!=minpoint1 && i!=minpoint2) 
                        { 
                            mindist = dist;  
                            minpoint3 = i; 
                        }
                    }
                if (mindist>maxdistall) { maxdistall = mindist; }

                // weighted average of 3 nearest points
                scalar w1 = 1/mag(p-coordinates[minpoint1]);
                scalar w2 = 1/mag(p-coordinates[minpoint2]);
                scalar w3 = 1/mag(p-coordinates[minpoint3]);
                // write nearest value to current patch
                interpHeat = (w1 * heat[minpoint1] 
                                + w2 * heat[minpoint2] 
                                + w3 * heat[minpoint3]) 
                                / (w1 + w2 + w3);
            } // method_
            QJ[cellsi] = interpHeat;

        } // cellsi
        forAll(QJ.boundaryField(), patchesi)
        {
            forAll(QJ.boundaryField()[patchesi], facesi)
            {
                vector p = mesh.Cf().boundaryField()[patchesi][facesi];
                vector por = p;
                // calculte r and project mesh points to 2D
                scalar rr = 0;
                // calculate radius rr from axis_
                if(rot2D_)
                {
                    forAll (p, xi)
                    {
                        rr = rr + (1-axis_[xi])*p[xi]*p[xi];
                    }
                    rr = Foam::sqrt(rr);
                    // projected p
                    p[rrLabel] = rr;
                    p[axisLabel] = p[axisLabel];
                    p[rotLabel] = 0;
                }

                scalar mindist = 1000;
                label minpoint = 0;

                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) { mindist = dist;  minpoint = i; }
                }

                if (mindist>maxdistall) { maxdistall = mindist; }
                
                scalar interpHeat(Zero);
                if (method_ =="nearest")
                {
                    interpHeat = heat[minpoint];
                }
                else if (method_ == "interp1")
                {
                    // find 2 nearest points to p
                    label minpoint1=0, minpoint2=0;
                    scalar mindist = 1000;
                    for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist) { mindist = dist;  minpoint1 = i; }
                    }
                    if (mindist>maxdistall) { maxdistall = mindist; }

                    mindist = 1000;
                    for (label i=0; i<coordinates.size(); i++) 
                        {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist && i!=minpoint1) { mindist = dist;  minpoint2 = i; }
                        }
                    if (mindist>maxdistall) { maxdistall = mindist; }

                    // weighted average of 2 nearest points
                    scalar dist1 = mag(p-coordinates[minpoint1]);
                    scalar dist2 = mag(p-coordinates[minpoint2]);
                    // write nearest value to current patch
                    scalar w = dist1/(dist1+dist2);
                    interpHeat = (1-w)*heat[minpoint1] + w*heat[minpoint2];
                }
                QJ.boundaryFieldRef()[patchesi][facesi] = interpHeat;
            }
        }
    }

    if(FLmap_)
    {
        DynamicList<vector> coordinates;
        DynamicList<vector> force;
        Info << "  Reading file "  << FLfile_ << endl;
        Info << "  coordinates at columns (x,y,z) = ["
            << coordColumn_ << ", " << coordColumn_+1 << ", " << coordColumn_+2 
            << "]" << endl;
        Info << "  Lorentz Force at columns (Fx,Fy) = ["
            << forceColumn_ << ", " << forceColumn_+1 << "]" 
            << endl;
            
        while (inFile.good())
        {
            string line;
            inFile.getLine(line);

            if (line.size()>1)
            {
                // Set up up a stream from this line
                std::stringstream lineStr(line);

                DynamicList<scalar> data;
                string item;
                // Info << "line = " << line << endl;

                std::string token;
                while (std::getline(lineStr, token, ',')) 
                {            
                    double value = std::stod(token);  // Convert the token to a double
                    data.append(value);      // Store the value in the vector
                }

                // read Elmer coordinates
                //TBD: add permutation
                vector coord(Zero);
                forAll (coord, xi)
                {
                    coord[xi] = data[coordColumn_+xi];
                }   

                // read elmer force
                vector forceVec(Zero);
                forceVec[0] = data[forceColumn_+0];
                forceVec[1] = data[forceColumn_+1];
                forceVec[2] = data[forceColumn_+2];
                
                coordinates.append(coord);
                force.append(forceVec);
            }
        }
        inFile.rewind();

        Info<< "  coordinates.size() = " << coordinates.size() << endl;
        Info<< "  force.size() = " << force.size() << endl;
        Info<< "  Interpolation method: " << method_ << endl;
        
        // Mapping
        Info<< "\nMapping data\n" << endl;

        scalar maxdistall = -1000;

        forAll(mesh.cells(), cellsi)
        {
            vector p = mesh.C()[cellsi];
            vector por = p;
            // calculte r
            scalar rr = 0;
            // calculate radius rr from axis_
            if(rot2D_)
            {
                forAll (p, xi)
                {
                    rr = rr + (1-axis_[xi])*p[xi]*p[xi];
                }
                rr = Foam::sqrt(rr);
                // flatten p
                p[rrLabel] = rr;
                p[axisLabel] = p[axisLabel];
                p[rotLabel] = 0;
            }

            // find nearest neighbour
            scalar mindist = 1000;
            label minpoint = 0;
            for (label i=0; i<coordinates.size(); i++) 
            {
                scalar dist = mag(p-coordinates[i]);
                if (dist < mindist) 
                { 
                    mindist = dist;  
                    minpoint = i; 
                }
            }
            // Info << "por = " << por << "; p = " << p << "; minpoint = " << coordinates[minpoint] << endl; 
            if (mindist>maxdistall) 
            { 
                maxdistall = mindist; 
            }     

            vector interpForce(Zero);

            if (method_ =="nearest")
            {
                interpForce = force[minpoint];
            }
            else if (method_ == "interp1")
            {
                // find 3 nearest points to p
                label minpoint1=0, minpoint2=0, minpoint3=0;
                // point1
                scalar mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) 
                    { 
                        mindist = dist;  
                        minpoint1 = i; 
                    }
                }
                if (mindist>maxdistall) { maxdistall = mindist; }
                // point2
                mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist && i!=minpoint1) 
                        { 
                            mindist = dist;  
                            minpoint2 = i; 
                        }
                    }
                if (mindist>maxdistall) { maxdistall = mindist; }

                // point3
                mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist && i!=minpoint1 && i!=minpoint2) 
                        { 
                            mindist = dist;  
                            minpoint3 = i; 
                        }
                    }
                if (mindist>maxdistall) { maxdistall = mindist; }

                // weighted average of 3 nearest points
                scalar w1 = 1/mag(p-coordinates[minpoint1]);
                scalar w2 = 1/mag(p-coordinates[minpoint2]);
                scalar w3 = 1/mag(p-coordinates[minpoint3]);
                // write nearest value to current patch
                interpForce = (w1 * force[minpoint1] 
                                + w2 * force[minpoint2] 
                                + w3 * force[minpoint3]) 
                                / (w1 + w2 + w3);
            } // method_

            vector vecRot(Zero);
            if(rot2D_)
            {
                if (p[rrLabel]>1e-10) 
                {
                    vecRot[rrLabel] = interpForce[rrLabel]*por[rrLabel]/(p[rrLabel]+SMALL);
                    vecRot[rotLabel] = interpForce[rrLabel]*por[rotLabel]/(p[rrLabel]+SMALL); 
                }
                else 
                {
                    vecRot[rrLabel] = 0;
                    vecRot[rotLabel] = 0;
                }
                
                vecRot[axisLabel] = interpForce[axisLabel];
            }
            else
            {
                vecRot = interpForce;
            }
            FL[cellsi] = vecRot;

        } // cellsi
        forAll(FL.boundaryField(), patchesi)
        {
            forAll(FL.boundaryField()[patchesi], facesi)
            {
                vector p = mesh.Cf().boundaryField()[patchesi][facesi];
                vector por = p;
                // calculte r
                scalar rr = 0;
                // calculate radius rr from axis_
                if(rot2D_)
                {
                    forAll (p, xi)
                    {
                        rr = rr + (1-axis_[xi])*p[xi]*p[xi];
                    }
                    rr = Foam::sqrt(rr);
                    // flatten p
                    p[rrLabel] = rr;
                    p[axisLabel] = p[axisLabel];
                    p[rotLabel] = 0;
                }

                // find nearest neighbour
                scalar mindist = 1000;
                label minpoint = 0;
                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) 
                    { 
                        mindist = dist;  
                        minpoint = i; 
                    }
                }

                vector interpForce(Zero);
                if (method_ =="nearest")
                {
                    interpForce = force[minpoint];
                }
                else if (method_ == "interp1")
                {
                    // find 3 nearest points to p
                    label minpoint1=0, minpoint2=0, minpoint3=0;
                    // point1
                    scalar mindist = 1000;
                    for (label i=0; i<coordinates.size(); i++) 
                    {
                        scalar dist = mag(p-coordinates[i]);
                        if (dist < mindist) 
                        { 
                            mindist = dist;  
                            minpoint1 = i; 
                        }
                    }
                    if (mindist>maxdistall) { maxdistall = mindist; }
                    // point2
                    mindist = 1000;
                    for (label i=0; i<coordinates.size(); i++) 
                        {
                            scalar dist = mag(p-coordinates[i]);
                            if (dist < mindist && i!=minpoint1) 
                            { 
                                mindist = dist;  
                                minpoint2 = i; 
                            }
                        }
                    if (mindist>maxdistall) { maxdistall = mindist; }

                    // point3
                    mindist = 1000;
                    for (label i=0; i<coordinates.size(); i++) 
                        {
                            scalar dist = mag(p-coordinates[i]);
                            if (dist < mindist && i!=minpoint1 && i!=minpoint2) 
                            { 
                                mindist = dist;  
                                minpoint3 = i; 
                            }
                        }
                    if (mindist>maxdistall) { maxdistall = mindist; }

                    // weighted average of 3 nearest points
                    scalar w1 = 1/mag(p-coordinates[minpoint1]);
                    scalar w2 = 1/mag(p-coordinates[minpoint2]);
                    scalar w3 = 1/mag(p-coordinates[minpoint3]);
                    // write nearest value to current patch
                    interpForce = (w1 * force[minpoint1] 
                                    + w2 * force[minpoint2] 
                                    + w3 * force[minpoint3]) 
                                    / (w1 + w2 + w3);
                } // method_

                vector vecRot(Zero);
                if(rot2D_)
                {
                    if (p[rrLabel]>1e-10) 
                    {
                        vecRot[rrLabel] = interpForce[rrLabel]*por[rrLabel]/(p[rrLabel]+SMALL);
                        vecRot[rotLabel] = interpForce[rrLabel]*por[rotLabel]/(p[rrLabel]+SMALL); 
                    }
                    else 
                    {
                        vecRot[rrLabel] = 0;
                        vecRot[rotLabel] = 0;
                    }
                    
                    vecRot[axisLabel] = interpForce[axisLabel];
                }
                else
                {
                    vecRot = interpForce;
                }
                FL.boundaryFieldRef()[patchesi][facesi] = vecRot;
            }
        }
    } // if FLmap_
} // if FLmap_ || QJmap_
