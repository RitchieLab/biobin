/*
 * Information.h
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_INFORMATION_H
#define KNOWLEDGE_INFORMATION_H

#include <string>
#include <ostream>
#include <vector>
#include <map>
#include <set>
#include <boost/iterator/iterator_facade.hpp>

using std::ostream;
using std::string;
using std::vector;
using std::map;
using std::set;

namespace Knowledge {

class Locus;
class Region;

/*!
 * \brief Defines a class that gets general information about the database.
 * This class is an interface to get specific information about the database.
 * The information returned by this class is meta-information about the
 * database, like the genomic build and the list of populations available.
 */
class Information {
public:

	/*!
	 * An enumeration that can be used as a bitmask to determine the role of a
	 * SNP.
	 */
	class snp_role {
	public:
		struct Ptr_Less: public std::binary_function<const snp_role*, const snp_role*, bool> {
			bool operator()(const snp_role* x, const snp_role* y) {
				return (y != 0 && x != 0) ? static_cast<int>(*x) < static_cast<int>(*y) : y < x;
			}
		};

		class const_iterator: public boost::iterator_facade<const_iterator,
				const snp_role&, boost::forward_traversal_tag> {

		public:
			const_iterator(set<const snp_role*, Ptr_Less>::const_iterator itr) : _itr(itr) {}

		private:
			friend class boost::iterator_core_access;

			void increment() {++_itr;}
			bool equal(const const_iterator& other) const {
				return _itr == other._itr;
			}
			const snp_role & dereference() const {return **_itr;}

			set<const snp_role*, Ptr_Less>::const_iterator _itr;
		};

	private:
		explicit snp_role(const string& val) :	_data(val) {
			if (s_val_map.find(_data) == s_val_map.end()) {
				s_val_map.insert(std::make_pair(_data, 1 << (++s_num_vals - 1)));
				s_enums.insert(this);
			}
		}



	public:

		friend class Information;

		operator int() const { return s_val_map[_data];}
		operator string() const {return _data;}

		static const_iterator begin() { return const_iterator(s_enums.begin());}
		static const_iterator end() { return const_iterator(s_enums.end()); }

	private:

		string _data;
		static int s_num_vals;
		static map<string, int> s_val_map;
		static set<const snp_role*, Ptr_Less> s_enums;
	};

	static const snp_role INTRON;
	static const snp_role EXON;
	static const snp_role REGULATORY;

	/*!
	 * Destroys the Information object
	 */
	virtual ~Information() {
	}

	/*!
	 * \brief Converts a population string into a Population ID.
	 * Takes the given pop_str and returns an integer ID needed for DB queries.
	 *
	 * \param pop_str The Population label.
	 *
	 * \return The integer index of the given population, or 0 if not found.
	 */
	virtual int getPopulationID(const string& pop_str) = 0;

	/*!
	 * \brief Gets the names of the sources for the given IDs
	 * Queries the database and returns a map of IDs to source name.
	 *
	 * \param group_ids A set of IDs to obtain sourse information for.  If empty,
	 * return source information for all available sources
	 * \param[out] type_names_out A mapping of ID->source name
	 */
	virtual void getGroupTypes(const set<unsigned int>& group_ids,
			map<int, string>& type_names_out) = 0;

	/*!
	 * \brief Gets the zone size.
	 * Queries the database and returns the size of the zones used in building
	 * the region_zone table.
	 */
	virtual int getZoneSize() = 0;

	/*!
	 * \brief Returns a SNP's role.
	 * This function returns a SNPs role in a gene.
	 *
	 * \param loc The Locus object in question
	 * \param reg The associated Region to get the role for.
	 *
	 * \return An integer that represents a bitmask of the snp_roles
	 */
	virtual int getSNPRole(const Locus& loc, const Region& reg) = 0;

	/*!
	 * \brief Prints a list of the populations.
	 * Prints a list of all of the populations available in the database
	 *
	 * \param os The output stream to print to.
	 */
	virtual void printPopulations(ostream& os) = 0;

	/*!
	 * \brief Prints a list of all of the available sources
	 * Prints all of the avaialable sources to the given output stream.
	 *
	 * \param os The output stream to print to.
	 */
	virtual void printSources(ostream& os) = 0;

	/*!
	 * \brief Returns a vector of the source IDs from the DB
	 * Returns the list of source IDs from the DB.  The source IDs are either
	 * those that correspond to the c_source_names file, or all source IDs in
	 * the database.
	 *
	 * \return A list of the numerical IDs of the sources
	 */
	virtual const set<unsigned int>& getSourceIds() = 0;

	/*!
	 * \brief Returns a string of IDs compatible with a "where" clause.
	 * Returns a string of IDs
	 */
	string getSourceList();

	static vector<string> c_source_names;

protected:
	static set<unsigned int> _s_source_ids;
};

}

#endif /* INFORMATION_H_ */
