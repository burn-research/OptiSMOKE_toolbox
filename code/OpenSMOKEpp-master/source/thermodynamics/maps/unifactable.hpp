/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

UnifacTable::UnifacTable() { }

void UnifacTable::SetupTable(const fs::path& map_path)
{
  lc::Config doc;
  try
    {
      int ngroups = 0;
      doc.readFile(map_path.string().c_str());
      lc::Setting &tableList = doc.lookup("unifac");
      lc::Setting &table = tableList["table"];
      ngroups = table.getLength();
      TableAllocation(ngroups);
      for (int i = 0; i < ngroups; i++)
        {
          table[i].lookupValue("subGroupName", groupname_[i]);
          table[i].lookupValue("mainGroup", maingroup_[i]);
          table[i].lookupValue("subGroup", subgroup_[i]);
          table[i].lookupValue("R", R_[i]);
          table[i].lookupValue("Q", Q_[i]);
        }
    }
	catch (const lc::SettingNotFoundException &nfex)
	{
		std::cerr	<< "Parse error at " << nfex.getPath() << std::endl;
	}
	catch (const lc::ParseException &pex)
	{
		std::cerr	<< "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
	}
	catch (const lc::FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file: " << fioex.what() << std::endl;
	}
}

void UnifacTable::TableAllocation(int ngroups)
{
  groupname_.resize(ngroups);
  maingroup_.resize(ngroups);
  subgroup_.resize(ngroups);
  R_.resize(ngroups);
  Q_.resize(ngroups);
}

void UnifacTable::SetupInteractionCoefficients(const fs::path& map_path)
{
  lc::Config doc;
  doc.readFile(map_path.string().c_str());

  try
    {
      lc::Setting &table = doc.lookup("InteractionCoefficients");
      int nmaingroups = table.getLength();
      a_.resize(nmaingroups, nmaingroups);
      for (int i = 0; i < nmaingroups; i++)
        {
          string name = "a" + OpenSMOKE::ToString(i+1);
          if (table[name.c_str()].getLength() != nmaingroups)
            {
              std::cerr << "Interaction matrix is not square. "
                      "Error at line " << i+1 << "/" << nmaingroups << std::endl;
              std::cerr << table[name.c_str()].getLength() <<
                      " elements" << std::endl;
              exit(-1);
            }
          for (int j = 0; j < nmaingroups; j++)
              a_(i, j) = table[name.c_str()][j];
        }
    }
  catch (const lc::SettingNotFoundException &nfex)
    {
      std::cerr << "Parse error at " << nfex.getPath() << std::endl;
    }
  catch (const lc::ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    }

}
