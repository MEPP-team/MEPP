///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_SCENE
#define HEADER_SCENE

#include <QObject>
#include <QFileInfo>

#include "./Polyhedron/polyhedron.h"

class Scene
{
	public:
		Scene();
		~Scene();

		int add_mesh(QString filename, int loadType, typeFuncOpenSave f, Viewer* viewer);
		int save_file(QString filename, typeFuncOpenSave f, Viewer* viewer);

		QString currentFile() { return curFile; }
		QString userFriendlyCurrentFile() { return strippedName(curFile); }

		void setcurrentFile(QString name) { curFile = name; }

		int get_nb_polyhedrons() { return (int)m_polyhedrons.size(); }

		int get_current_polyhedron() { return m_current_polyhedron; }
		void set_current_polyhedron(int current_polyhedron) { m_current_polyhedron = current_polyhedron; }

		PolyhedronPtr const get_polyhedron() { return m_polyhedrons.at(get_current_polyhedron()); }
		PolyhedronPtr const get_polyhedron(unsigned int i) { return m_polyhedrons.at(i); }

		void add_polyhedron(PolyhedronPtr polyhedron_ptr) { m_polyhedrons.push_back(polyhedron_ptr); }
		void delete_polyhedron(unsigned int i) { m_polyhedrons.erase(m_polyhedrons.begin()+i); }

		vector<PolyhedronPtr>::iterator get_begin_polyhedrons() { return m_polyhedrons.begin(); }
		void insert_polyhedron(vector<PolyhedronPtr>::iterator it, PolyhedronPtr polyhedron_ptr) { m_polyhedrons.insert(it, polyhedron_ptr); }

		int get_loadType() { return m_loadType; }
		void set_loadType(int loadType) { m_loadType = loadType; }
		QString get_stringLoadType()
		{
			switch (m_loadType)
			{
				case 0: return "Normal";
				case 1: return "Space";
				case 2: return "Time";
				default: return "Normal";
			}
		}

	private:
		vector<PolyhedronPtr> m_polyhedrons;
		int m_current_polyhedron;
		int m_loadType;

		int load_file(PolyhedronPtr polyhedron_ptr, QString filename);

		//

		QString curFile;
		QString strippedName(const QString &fullFileName) { return QFileInfo(fullFileName).fileName(); }
};

typedef boost::shared_ptr<Scene> ScenePtr;

#endif
