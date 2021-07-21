%> @file  mainPioneer.m
%> @brief Documentation for image reconstruction Pioneer
%======================================================================
%> @mainpage Documentation for image reconstruction Pioneer
%>
%> @section intro Introduction
%>
%> The @b Image @b Reconstruction @b Pioneer is a package based on MATLAB for 
%> near infrared optical tomography.
%>
%> It incorporates various image reconstruction methods, forward models and
%> tissue geometries, source detector arrangements, etc.
%>
%> @section Forward models
%> Numerical models include:
%>
%> - Monte Carlo Methods
%> - Finite Element Methods
%>
%> 
%>   @note for enumeration definition : this script only supports the following declaration :
%>   @code
%>   enumeration
%>     first_enum
%>     second_enum
%>   end
%>   @endcode
%>   This one is not yet supported :
%>   @code
%>   enumeration
%>     first_enum, second_enum
%>   end
%>   @endcode
%>   The file classDocumentationExample.m provides an example of enumeration definition
%>   with comments extracted by Doxygen.
%>
%> .
%>
%> See README.md if you want more details about how to make this script work.
%>
%> @attention Each line belonging to the doxygen documentation must begin with <b>%></b> .
%>
%> @par Example
%>
%>@code
%>% Matlab comment ignored by doxygen
%>%> comment analyzed by doxygen
%>@endcode
%> @attention Doxygen keyword have to begin by @b @, for example @@b to bold the text (the use of \ instead of @ is not supported)
%>
%>@section funcDecr Function description
%> The keyword @b @@param and @b @@retval will be used to describe the input and
%> output parameters of a function.
%>
%> For function description, the description should follow the following presentation :
%>
%>@verbatim
%>% ======================================================================
%>%> @brief Brief description of the function
%>%>
%>%> More detailed description.
%>%>
%>%> @param arg1 First argument
%>%> @param arg2 Second argument
%>%>
%>%> @retval out1 return value for the first output variable
%>%> @retval out2 return value for the second output variable
%>% ======================================================================
%>[out1, out2] = function( arg1, arg2)
%>  out1 = arg2;
%>  out2 = arg1;
%>end
%>@endverbatim
%>
%>
%> @section classDecr Class description
%>
%> For class description, the following description can be used :
%>
%>@verbatim
%>% ======================================================================
%>%> @brief Brief description of the class
%>%>
%>%> More detailed description of the class
%>%>
%>% ======================================================================
%>@endverbatim
%>
%>

%======================================================================
%> @brief Brief description of the function
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================

%%