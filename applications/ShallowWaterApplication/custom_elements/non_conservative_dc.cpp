//
//   Project Name:        Kratos
//   Last modified by:    $Author:  Miguel Masó Sotomayor $
//   Date:                $Date:             june 14 2017 $
//   Revision:            $Revision:                  1.1 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/non_conservative_dc.h"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	NonConservativeDC::NonConservativeDC(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		// DO NOT ADD DOFS HERE
	}

	//************************************************************************************
	//************************************************************************************
	NonConservativeDC::NonConservativeDC(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer NonConservativeDC::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new NonConservativeDC(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	NonConservativeDC::~NonConservativeDC()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void NonConservativeDC::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
	{
		const unsigned int number_of_nodes = 3;
		//~ const unsigned int number_of_dof = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		for(unsigned int i = 0; i<number_of_nodes; i++){
			for(unsigned int j = 0; j<number_of_nodes; j++){
				for(unsigned int k = 0; k<3; k++)
					rMassMatrix(k+3*i, k+3*j) += 1.0;
			}
		}
		rMassMatrix *= 1.0/(12.0);
	}

	void NonConservativeDC::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
	{
		const unsigned int number_of_nodes = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		rMassMatrix *= 1.0/(3.0);
	}

	//************************************************************************************
	//************************************************************************************
	void NonConservativeDC::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		//~ const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		//~ array_1d<double,2> BDFcoeffs = {1.0, 1.0};  
		double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};

		boost::numeric::ublas::bounded_matrix<double,9,9> msMass   = ZeroMatrix(9,9);     // Mass matrix
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX  = ZeroMatrix(3,2);     // Shape functions gradients
		boost::numeric::ublas::bounded_matrix<double,9,9> msC      = ZeroMatrix(9,9);     // Nt*A*B (LHS)
		boost::numeric::ublas::bounded_matrix<double,3,9> msN      = ZeroMatrix(3,9);     // Shape functions type
		//
		boost::numeric::ublas::bounded_matrix<double,2,9> msN_vel        = ZeroMatrix(2,9);   // Shape functions type matrix (for velocity unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msN_height     = ZeroMatrix(1,9);   // Shape functions type matrix (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,2,9> msDN_DX_height = ZeroMatrix(2,9);   // Shape functions gradients matrix (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msDN_DX_vel    = ZeroMatrix(1,9);   // Shape functions gradients matrix (for velocity unknown)
		//
		array_1d<double,3> msNGauss;                                    // Dimension = number of nodes . Position of the gauss point
		array_1d<double,9> ms_depth;
		array_1d<double,9> ms_unknown;
		array_1d<double,9> ms_proj_unknown;
		double height;
		double k_dc;                                                                       // Discontinuity capturing

		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points*3)
			rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false); // Resizing the system in case it does not have the right size
		if(rRightHandSideVector.size() != number_of_points*3)
			rRightHandSideVector.resize(number_of_points*3,false);

		// Getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msNGauss, Area); // Asking for gradients and other info
		double elem_size = pow(Area,0.5);

		// Reading properties and conditions
		//~ const double gravity = rCurrentProcessInfo[GRAVITY_Z];
		const double gravity = 9.81;

		// Get current step and projected values
		int counter = 0;
		for(unsigned int iii = 0; iii<number_of_points; iii++){
			ms_depth[counter] = 0;
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_X);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
			counter++;

			ms_depth[counter] = 0;
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_Y);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
			counter++;

			ms_depth[counter]   = GetGeometry()[iii].FastGetSolutionStepValue(BATHYMETRY);
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_HEIGHT);
			counter++;
		}

		// Compute parameters and derivatives matrices
		// Loop on Gauss points: ONE GAUSS POINT

		// Height gradient
		msDN_DX_height(0,2) = msDN_DX(0,0);
		msDN_DX_height(0,5) = msDN_DX(1,0);
		msDN_DX_height(0,8) = msDN_DX(2,0);
		msDN_DX_height(1,2) = msDN_DX(0,1);
		msDN_DX_height(1,5) = msDN_DX(1,1);
		msDN_DX_height(1,8) = msDN_DX(2,1);
		// Velocity divergence
		msDN_DX_vel(0,0) = msDN_DX(0,0);
		msDN_DX_vel(0,1) = msDN_DX(0,1);
		msDN_DX_vel(0,3) = msDN_DX(1,0);
		msDN_DX_vel(0,4) = msDN_DX(1,1);
		msDN_DX_vel(0,6) = msDN_DX(2,0);
		msDN_DX_vel(0,7) = msDN_DX(2,1);
		
		// N matrix: shape functions
		for (unsigned int jjj = 0; jjj < number_of_points; jjj++) {
			for (unsigned int iii = 0; iii < number_of_points; iii++)
				msN(iii,iii+3*jjj) = msNGauss[jjj];
		}
		//noalias(msNmass) = row(msN, 2);
		msN_height(0,2) = msNGauss[0];
		msN_height(0,5) = msNGauss[1];
		msN_height(0,8) = msNGauss[2];
		//noalias(msNmoment) = row(msN, 0,1);
		msN_vel(0,0) = msNGauss[0];
		msN_vel(0,3) = msNGauss[1];
		msN_vel(0,6) = msNGauss[2];
		msN_vel(1,1) = msNGauss[0];
		msN_vel(1,5) = msNGauss[1];
		msN_vel(1,7) = msNGauss[2];

		height = msNGauss[0]*ms_unknown[2] + msNGauss[1]*ms_unknown[5] + msNGauss[2]*ms_unknown[8];

		// End loop on Gauss point


		// Main loop
		// LHS
		// Cross terms
		noalias(rLeftHandSideMatrix)  = height*prod(trans(msN_height),msDN_DX_vel);        // Mass: q*h*div(u)
		noalias(msC)                  = gravity*prod(trans(msN_vel),msDN_DX_height);       // Momentum: w*g*grad(h)
		noalias(rLeftHandSideMatrix) += msC;

		// Inertia terms
		// Compute the mass matrix
		CalculateConsistentMassMatrix(msMass);
		//~ CalculateLumpedMassMatrix(msMass);
		// LHS += bdf*M
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMass;

		// Add Discontinuity capturing term via adding artificial diffusion to velocity
		//~ double divergence_is_zero = 1e-1;
		//~ double max_viscosity = 100;
		//~ double m_residual;
		//~ double m_divergence;
		//~ m_residual = norm_2(gravity*prod(msDN_DX_height, ms_depth - ms_unknown) + BDFcoeffs[1]*prod(msN_vel, (ms_proj_unknown - ms_unknown) ) );
		//~ this->SetValue(RESIDUAL_NORM,m_residual);
		//~ m_divergence = norm_2(prod(msDN_DX_vel, ms_unknown));
		//~ if (m_divergence < divergence_is_zero){
			//~ k_dc = 0;
			//~ m_divergence = 0;
		//~ }
		//~ else {
			//~ k_dc = 0.1*0.5*m_residual/m_divergence;
		//~ }
		//~ this->SetValue(MIU,m_divergence);
		//~ if (k_dc > max_viscosity){
			//~ k_dc = 0;
		//~ }
		//~ this->SetValue(VEL_ART_VISC,k_dc);
		//~ noalias(rLeftHandSideMatrix) += k_dc * prod(trans(msDN_DX_vel), msDN_DX_vel);

		// Add discontinuity capturing term via adding artificial diffusion to height
		double gradient_is_zero = 1e-6;
		
		// A. Artificial diffusion via RESIDUAL
		double m_residual;
		double m_height_grad_norm;
		//~ 
		m_residual = norm_1(prod(msN_height,ms_unknown)) * norm_1(prod(msDN_DX_vel,ms_unknown)) + BDFcoeffs[1]*norm_1(prod(msN_height, (ms_unknown - ms_proj_unknown)));
		m_height_grad_norm = norm_2(prod(msDN_DX_height,ms_unknown));
		if (m_height_grad_norm < gradient_is_zero){
			k_dc = 0;
		}
		else{
			k_dc = 0.5*0.4*elem_size*m_residual;//m_height_grad_norm;  // Residual formulation
		}
		
		
		// B. Artificial diffusion via "VELOCITY PROJECTION"
		//~ array_1d<double,2> m_velocity;
		//~ array_1d<double,2> m_height_grad;
		//~ double m_height_grad_norm;
		//~ double m_vel_project_norm;
		//~ 
		//~ m_velocity = prod(msN_vel,ms_unknown);
		//~ m_height_grad = prod(msDN_DX_height,ms_unknown);
		//~ m_height_grad_norm = norm_2(m_height_grad);
		//~ m_vel_project_norm = abs(inner_prod(m_velocity,m_height_grad));
		//~ if (m_height_grad_norm < gradient_is_zero){
			//~ k_dc = 0;
		//~ }
		//~ else{
			//~ k_dc = 0.5*0.1*elem_size*m_vel_project_norm/m_height_grad_norm;
		//~ }
		
		
		// C. Artificial diffusion via ORTHOGONAL PROJECTION
		//~ array_1d<double,2> m_velocity;
		//~ array_1d<double,2> m_height_grad;
		//~ double m_height_grad_norm;
		//~ double m_vel_norm;
		//~ double m_orthogonalP;
		//~ boost::numeric::ublas::bounded_matrix<double,2,9> mP_height;
		//~ 
		//~ m_velocity    = prod(msN_vel,ms_unknown);
		//~ m_height_grad = prod(msDN_DX_height,ms_unknown);
		//~ m_vel_norm         = norm_2(m_velocity);
		//~ m_height_grad_norm = norm_2(m_height_grad);
		//~ 
		//~ for (unsigned int iii = 0; iii<3*number_of_points; iii++){
			//~ mP_height(0,iii) = msDN_DX_height(0,iii) * msDN_DX_height(0,iii); 
			//~ mP_height(1,iii) = msDN_DX_height(1,iii) * msDN_DX_height(1,iii);
		//~ }
		//~ m_orthogonalP = norm_2(m_height_grad - prod(mP_height,ms_unknown));
		//~ 
		//~ if (m_height_grad_norm < gradient_is_zero){
			//~ k_dc = 0;
		//~ }
		//~ else{
			//~ k_dc = 0.5*0.2*elem_size*m_vel_norm*m_orthogonalP/m_height_grad_norm;
		//~ }
		//~ this->SetValue(RESIDUAL_NORM,m_orthogonalP);
		//~ this->SetValue(MIU,m_height_grad_norm);
		//~ this->SetValue(PR_ART_VISC,k_dc);
		//~ KRATOS_WATCH(m_vel_project_norm)
		//~ KRATOS_WATCH(m_height_grad_norm)
		//~ KRATOS_WATCH(k_dc)
		//~ KRATOS_WATCH(k_dc * prod(trans(msDN_DX_height), msDN_DX_height))
		
		
		// D. My own artificial diffusion
		//~ array_1d<double,2> m_height_grad;
		//~ double m_height_grad_norm;
		//~ m_height_grad = prod(msDN_DX_height,ms_unknown);
		//~ m_height_grad_norm = norm_2(m_height_grad);
		//~ if (m_height_grad_norm < gradient_is_zero){
			//~ k_dc = 0;
		//~ }
		//~ else{
			//~ k_dc = 0.5*0.4*elem_size*m_height_grad_norm;
		//~ }
		this->SetValue(MIU,m_height_grad_norm);
		this->SetValue(RESIDUAL_NORM,m_residual);
		this->SetValue(PR_ART_VISC,k_dc);
		
		noalias(rLeftHandSideMatrix) += k_dc * prod(trans(msDN_DX_height), msDN_DX_height);

		// RHS
		// TODO: SOURCE TERM
		noalias(rRightHandSideVector) = -prod(msC, ms_depth);

		// Inertia terms
		// RHS += M*vhistory
		noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(msMass, ms_proj_unknown);

		// Substracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_unknown);

		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void NonConservativeDC::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	// This subroutine calculates the nodal contributions for the explicit steps of the
	// Fractional step procedure
	void NonConservativeDC::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void NonConservativeDC::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes*3)
			rResult.resize(number_of_nodes*3,false);
		int counter=0;
		
		for (unsigned int i=0;i<number_of_nodes;i++){
			rResult[counter++] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(HEIGHT).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void NonConservativeDC::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rElementalDofList.size() != number_of_nodes*3)
			rElementalDofList.resize(number_of_nodes*3);
		
		int counter=0;
		
		for (unsigned int i=0;i<number_of_nodes;i++){
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(VELOCITY_X);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(VELOCITY_Y);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(HEIGHT);
		}
	}

	void NonConservativeDC::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == VEL_ART_VISC){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
				rValues[PointNumber] = double(this->GetValue(VEL_ART_VISC));
		}
		if (rVariable == PR_ART_VISC){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
				rValues[PointNumber] = double(this->GetValue(PR_ART_VISC));
		}
		if (rVariable == RESIDUAL_NORM){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
				rValues[PointNumber] = double(this->GetValue(RESIDUAL_NORM));
		}
		if (rVariable == MIU){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
				rValues[PointNumber] = double(this->GetValue(MIU));
		}
    }
	//~ {
		//~ if (rVariable == MU)
		//~ {
			//~ // Shape functions and integration points
			//~ ShapeFunctionDerivativesArrayType DN_DX;
			//~ Matrix NContainer;
			//~ VectorType GaussWeights;
			//~ this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
			//~ const unsigned int NumGauss = GaussWeights.size();
//~ 
			//~ rValues.resize(NumGauss);
//~ 
			//~ double Density = 0.0;
//~ 
			//~ // Loop on integration points
			//~ for (unsigned int g = 0; g < NumGauss; g++)
			//~ {
				//~ const ShapeFunctionsType& N = row(NContainer,g);
				//~ const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
//~ 
				//~ this->EvaluateInPoint(Density,DENSITY,N);
				//~ double ElemSize = this->ElementSize();
//~ 
				//~ rValues[g] = this->EffectiveViscosity(Density,N,rDN_DX,ElemSize,rCurrentProcessInfo);
			//~ }
//~ 
		//~ }
		//~ else
		//~ {
		//~ }
	//~ }



} // namespace Kratos
