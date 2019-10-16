/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Hypre_helper

Group

Description

SourceFiles
    Hypre_helper.C

\*---------------------------------------------------------------------------*/
#ifndef HYPRE_HELPER_H
#define HYPRE_HELPER_H

#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//namespace Foam
//{

/*---------------------------------------------------------------------------*\
                             Class Hypre_helper Declaration
\*---------------------------------------------------------------------------*/

class Hypre_helper
{

public:

        //- Disallow default bitwise copy construct
        //Hypre_helper(const Hypre_helper&);

        //- Disallow default bitwise assignment
        //void operator=(const Hypre_helper&);

	bool first;

	// Hypre handles
	HYPRE_Solver solver;//, precond;
	HYPRE_IJMatrix A;
	HYPRE_ParCSRMatrix parcsr_A;
        int ilower, iupper, jlower, jupper;
	

    // Constructors

        //- Construct from matrix components and solver controls
        Hypre_helper():first(true){};


	//- Destructor
	~Hypre_helper()
	{
	   if(!first)
	   {
	   	HYPRE_BoomerAMGDestroy(solver);
	   	HYPRE_IJMatrixDestroy(A);
	   }
	};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
