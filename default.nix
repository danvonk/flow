{ lib
, llvmPackages_11
, cmake
, spdlog
, abseil-cpp
, cudaPackages
, vtk
, boost }:

cudaPackages.backendStdenv.mkDerivation rec {
  pname = "tflow";
  version = "0.1.0";
  
  src = ./.;

  nativeBuildInputs = [ cmake
        cudaPackages.cudatoolkit
                      ];
  buildInputs = [ spdlog
                  abseil-cpp
                  cudaPackages.cudatoolkit
                  vtk
                  boost


                ];

  cmakeFlags = [
    "-DENABLE_TESTING=OFF"
    "-DENABLE_INSTALL=ON"
  ];

  meta = with lib; {
    homepage = "https://github.com/nixvital/nix-based-cpp-starterkit";
    description = ''
      A template for Nix based C++ project setup.";
    '';
    licencse = licenses.mit;
    platforms = with platforms; linux ++ darwin;
    maintainers = [ maintainers.breakds ];    
  };
}
