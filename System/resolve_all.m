function [fList, pList] = resolve_all(Name)

  %find all functions needed by a routine
  %feed routine name to function as a string, returns a cell array of required files
  %
  %just remaps an internal Matlab capability to an easier-to-remember function name
  %Corwin Wright, c.wright@bath.ac.uk, 15/MAR/2016

  [fList, pList] = matlab.codetools.requiredFilesAndProducts(Name);

end

