/**
TODO: Rename this! 
    Why? Because it is confusing since it only generates local orbitlists from a given primitive orbitlist and a supercell
    It sounds like it can generate any kind of orbitlist


This is a small class or struct that will have a:

orbitlist (from primitive structure)
supercell
list of unique primitive cell offsets that the supercell span
the primToSupercellMap


you can query this object with

///Give number of possible local orbitlists one can make from current supercell and the primitive orbitlist
size_t number_of_possible_local_orbitlists();

///Generate the orbitlist from the primitive offset with count i
OrbitList getLocalOrbitlist(int i);

std::vector<Vector3d> getUniqueOffsets() const;

std::vector<Vector3d> primToSupercellMap() const;

///clears primToSupercellMap and unique offsets
void reset();

etc...
*/