/*!
 * \file scene.h
 * \brief Scene file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */ 
#ifndef HEADER_SCENE
#define HEADER_SCENE

#include <QObject>
#include <QFileInfo>

#include "./Polyhedron/polyhedron.h"

/*! 
 * \class Scene
 * \brief Scene class.
 */
class Scene
{
	public:
		/*!
		 * \brief Constructor.
		 */
		Scene();
		/*!
		 * \brief Destructor.
		 */
		~Scene();

		/*!
		 * \fn int add_mesh(QString filename, int loadType, typeFuncOpenSave f, Viewer* viewer)
		 * \brief Open and load a mesh in a viewer.
		 *
		 * \param filename the mesh filename.
		 * \param loadType load mode (Normal, Space or Time).
		 * \param f a pointer to a typeFuncOpenSave function to load the file (default function if NULL).
		 * \param viewer the viewer.
		 * \return 0 if loading is ok.
		 */
		int add_mesh(QString filename, int loadType, typeFuncOpenSave f, Viewer* viewer);
#ifdef WITH_ASSIMP
		/*!
		 * \fn int add_meshes_with_assimp(QString filename, Viewer* viewer)
		 * \brief Open and load meshes in a viewer with assimp library (if the file contain more than one mesh the load mode is set to Space)
		 *
		 * \param filename the mesh filename.
		 * \param viewer the viewer.
		 * \return 0 if loading is ok.
		 */
		int add_meshes_with_assimp(QString filename, Viewer* viewer);
#endif
		/*!
		 * \fn todoIfModeSpace(Viewer* viewer, double ystep)
		 * \brief Do everything we need for the space mode after a mesh loading (adding frame...).
		 *
		 * \param viewer the viewer.
		 * \param ystep 'y step' between meshes.
		 */
		void todoIfModeSpace(Viewer* viewer, double ystep);
		/*!
		 * \fn int save_file(QString filename, typeFuncOpenSave f, Viewer* viewer)
		 * \brief Save polyhedron(s) from the viewer in file(s).
		 *
		 * \param filename the mesh filename.
		 * \param f a pointer to a typeFuncOpenSave function to save the file(s) (default function if NULL).
		 * \param viewer the viewer.
		 * \return 0 if loading is ok.
		 */
		int save_file(QString filename, typeFuncOpenSave f, Viewer* viewer);

		/*!
		 * \fn QString currentFile()
		 * \brief Return filename of the current active polyhedron.
		 *
		 * \return QString.
		 */
		QString currentFile() { return QString::fromStdString(get_polyhedron()->pName); }
		/*!
		 * \fn QString userFriendlyCurrentFile()
		 * \brief Return filename (without path) of the current active polyhedron.
		 *
		 * \return QString.
		 */
		QString userFriendlyCurrentFile() { return strippedName(currentFile()); }
		/*!
		 * \fn QString currentFile(int p)
		 * \brief Return filename of the polyhedron of indice p.
		 *
		 * \param p polyhedron indice.
		 * \return QString.
		 */
		QString currentFile(int p) { return QString::fromStdString(get_polyhedron(p)->pName); }
		/*!
		 * \fn QString userFriendlyCurrentFile(int p)
		 * \brief Return filename (without path) of the polyhedron of indice p.
		 *
		 * \param p polyhedron indice.
		 * \return QString.
		 */
		QString userFriendlyCurrentFile(int p) { return strippedName(currentFile(p)); }

		/*!
		 * \fn setcurrentFile(QString name)
		 * \brief Set a new name for the current active polyhedron.
		 *
		 * \param name new polyhedron name.
		 */
		void setcurrentFile(QString name) { get_polyhedron()->pName = name.toStdString(); }
		/*!
		 * \fn setVisible(bool visible)
		 * \brief Toggle current active polyhedron visibility.
		 *
		 * \param visible true or false.
		 */
		void setVisible(bool visible) { get_polyhedron()->pShow = visible; }
		/*!
		 * \fn setcurrentFile(int p, QString name)
		 * \brief Set a new name for the polyhedron of indice p.
		 *
		 * \param p polyhedron indice.
		 * \param name new polyhedron name.
		 */
		void setcurrentFile(int p, QString name) { get_polyhedron(p)->pName = name.toStdString(); }
		/*!
		 * \fn setVisible(int p, bool visible)
		 * \brief Toggle visibility for polyhedron of indice p.
		 *
		 * \param p polyhedron indice.
		 * \param visible true or false.
		 */
		void setVisible(int p, bool visible) { get_polyhedron(p)->pShow = visible; }

		/*!
		 * \fn get_nb_polyhedrons()
		 * \brief Return the number of polyhedron(s) in the scene.
		 *
		 * \return int.
		 */
		int get_nb_polyhedrons() { return (int)m_polyhedrons.size(); }

		/*!
		 * \fn get_current_polyhedron()
		 * \brief Return the indice of the current active polyhedron in the scene.
		 *
		 * \return int.
		 */
		int get_current_polyhedron() { return m_current_polyhedron; }
		/*!
		 * \fn set_current_polyhedron(int current_polyhedron)
		 * \brief Set the current active polyhedron in the scene.
		 *
		 * \param current_polyhedron polyhedron indice.
		 */
		void set_current_polyhedron(int current_polyhedron) { m_current_polyhedron = current_polyhedron; }

		/*!
		 * \fn get_polyhedron()
		 * \brief Return the pointer of the current active polyhedron in the scene.
		 *
		 * \return PolyhedronPtr.
		 */
		PolyhedronPtr const get_polyhedron() { return m_polyhedrons.at(get_current_polyhedron()); }
		/*!
		 * \fn get_polyhedron(unsigned int i)
		 * \brief Return the pointer of the polyhedron of indice i in the scene.
		 *
		 * \param i polyhedron indice.
		 * \return PolyhedronPtr.
		 */
		PolyhedronPtr const get_polyhedron(unsigned int i) { return m_polyhedrons.at(i); }

		/*!
		 * \fn add_polyhedron(PolyhedronPtr polyhedron_ptr)
		 * \brief Add a polyhedron to the scene.
		 *
		 * \param polyhedron_ptr polyhedron pointeur.
		 */
		void add_polyhedron(PolyhedronPtr polyhedron_ptr) { m_polyhedrons.push_back(polyhedron_ptr); }
		/*!
		 * \fn delete_polyhedron(unsigned int i)
		 * \brief Delete a polyhedron from the scene.
		 *
		 * \param i polyhedron indice.
		 */
		void delete_polyhedron(unsigned int i) { m_polyhedrons.erase(m_polyhedrons.begin()+i); }

		/*!
		 * \fn get_begin_polyhedrons()
		 * \brief Return the beginning iterator of the list of polyhedrons.
		 *
		 * \return beginning iterator.
		 */
		vector<PolyhedronPtr>::iterator get_begin_polyhedrons() { return m_polyhedrons.begin(); }
		/*!
		 * \fn insert_polyhedron(vector<PolyhedronPtr>::iterator it, PolyhedronPtr polyhedron_ptr)
		 * \brief Insert a polyhedron to the scene.
		 *
		 * \param polyhedron_ptr polyhedron pointeur.
		 * \param it polyhedron iterator where to insert polyhedron just after.
		 */
		void insert_polyhedron(vector<PolyhedronPtr>::iterator it, PolyhedronPtr polyhedron_ptr) { m_polyhedrons.insert(it, polyhedron_ptr); }

		/*!
		 * \fn get_loadType()
		 * \brief Return the scene mode (Normal, Space or Time).
		 *
		 * \return int.
		 */
		int get_loadType() { return m_loadType; }
		/*!
		 * \fn set_loadType(int loadType)
		 * \brief Set the scene mode (Normal, Space or Time).
		 *
		 * \param loadType scene mode (Normal, Space or Time).
		 */
		void set_loadType(int loadType) { m_loadType = loadType; }
		/*!
		 * \fn get_stringLoadType()
		 * \brief Return a string of the scene mode ("Normal", "Space" or "Time").
		 *
		 * \return QString.
		 */
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
		/*!
		 * \fn get_stringLoadTypeForSwap()
		 * \brief Return a string of the available switching scene mode ("Normal", "-> to Time" or "-> to Space").
		 *
		 * \return QString.
		 */
		QString get_stringLoadTypeForSwap()
		{
			switch (m_loadType)
			{
				case 0: return "Normal";
				case 1: return "-> to Time";
				case 2: return "-> to Space";
				default: return "Normal";
			}
		}

	private:
		vector<PolyhedronPtr> m_polyhedrons;	//!< list of polyhedrons
		int m_current_polyhedron;				//!< indice of current active polyhedron
		int m_loadType;							//!< scene mode (Normal, Space or Time)

		//int load_file(PolyhedronPtr polyhedron_ptr, QString filename);

		//

		/*!
		 * \fn strippedName(const QString &fullFileName)
		 * \brief Return the filename of a mesh without path.
		 *
		 * \param fullFileName filename of the mesh.
		 * \return QString.
		 */
		QString strippedName(const QString &fullFileName) { return QFileInfo(fullFileName).fileName(); }
};

typedef boost::shared_ptr<Scene> ScenePtr;	//!< boost Scene shared pointeur

#endif
