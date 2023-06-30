# Provides access to the shell defined in flake.nix
(builtins.getFlake ("git+file://" + toString ./.)).devShells.${builtins.currentSystem}.default