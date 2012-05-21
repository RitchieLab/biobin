/*
 * RegionCollectionSQLite.h
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_REGIONCOLLECTIONSQLITE_H
#define KNOWLEDGE_REGIONCOLLECTIONSQLITE_H_

#include "RegionCollection.h"
#include "InformationSQLite.h"

#include <sqlite3.h>
#include "any_iterator.hpp"

namespace Knowledge{

/*!
 * \brief A class that implements the region collection from a SQLite database.
 * This class represents a RegionCollection from an SQLite LOKI database.
 * Users may construct this from either an already open SQLite3 database
 * connection, or from a filename that gives the location of the SQLite3 db.
 * the only dependency required for this class is the sqlite3 C interface
 * library (v. >=3.5.6), which should come standard with most modern Linux
 * distributions.
 */
class RegionCollectionSQLite : public RegionCollection{

	class Container {

	public:
		typedef IteratorTypeErasure::any_iterator<Knowledge::Locus* const, boost::bidirectional_traversal_tag> const_iterator;

		virtual ~Container() {}
		virtual const_iterator begin() const = 0;
		virtual const_iterator end() const = 0;
	};

	template <class T_cont>
	class LocusContainer : virtual public Container {
	public:

		LocusContainer(const T_cont& c) : _data(c) {}
		virtual ~LocusContainer() {}

		virtual const_iterator begin() const {return static_cast<const_iterator>(_data.begin());}
		virtual const_iterator end() const {return static_cast<const_iterator>(_data.end());}

	private:
		const T_cont& _data;
	};


public:
	/*!
	 * \brief Create a RegionCollection giving the DB location.
	 * Creates a RegionCollection from the fiven location of the sqlite3
	 * database.  This file must exist and conform to the LOKI specifications.
	 *
	 * \param fn The filename of the LOKI database.
	 */
	template <class T_cont>
	RegionCollectionSQLite(const string& fn, const T_cont& loci);
	/*!
	 * \brief Create a RegionCollection giving an open DB connection.
	 * Creates a RegionCollection with an already open databse connection that
	 * points to a LOKI database.  In this case, the user is responsible for
	 * closing the connection, and the connection cannot be closed until this
	 * object is destroyed.
	 *
	 * \param db The database connection to the LOKI database.
	 */
	template <class T_cont>
	RegionCollectionSQLite(sqlite3 *db, const T_cont& loci);
	/*!
	 * Destroys the RegionCollection objects, and if necessary closes the
	 * database connection.
	 */
	virtual ~RegionCollectionSQLite();

	/*!
	 * \brief Loading function implementation.
	 * This function loads all of the Regions that are identified by the ids or
	 * by the aliasList from the LOKI database into memory.  If both ids and
	 * aliasList are empty, will load all Regions from the databse into memory.
	 *
	 * \param ids A list of IDs to filter the results by.
	 * \param aliasList A list of aliases to filter results
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	virtual uint Load(const unordered_set<uint>& ids,
			const vector<string>& aliasList);

private:
	//! true if we opened the connection, false otherwise
	bool self_open;
	//! sqlite connection
	sqlite3 *db;

	//! object to get generalized information
	Information* _info;

	const Container* const _dataset;

	//! Adds a region based on the row (or returns the already added region)
	Knowledge::Region* addRegion(sqlite3_stmt* row);

	/*!
	 * Callback to parse a list of region IDs and add them to the undordered_list
	 */
	static int parseRegionIDQuery(void*, int, char**, char**);
};

template<class T_cont>
RegionCollectionSQLite::RegionCollectionSQLite(const string& fn,
		const T_cont& loci) :
		self_open(true), _dataset(new LocusContainer<T_cont>(loci)) {

	_info = new InformationSQLite(db);
	sqlite3_open(fn.c_str(), &db);
}

template<class T_cont>
RegionCollectionSQLite::RegionCollectionSQLite(sqlite3* db_conn,
		const T_cont& loci) :
		self_open(false), db(db_conn), _dataset(new LocusContainer<T_cont>(loci)) {
	_info = new InformationSQLite(db);
}

}




#endif /* REGIONCOLLECTIONSQLITE_H_ */
