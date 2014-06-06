/*!
 * \file viewer.hxx
 * \brief Viewer file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */  
#ifndef HEADER_VIEWER
#define HEADER_VIEWER

//#include <GL/glew.h>
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/manipulatedFrame.h> // fix for QGLViewer 2.5.2

//#include <QtGui/QCloseEvent>
#if QT_VERSION >= 0x050000
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#else
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#endif

#include "scene.h"

#ifdef WITH_FFMPEG
#include "QTFFmpegWrapper/QVideoEncoder.h"
#endif

/*! 
 * \class MeppManipulatedFrame
 * \brief MeppManipulatedFrame class.
 */
class MeppManipulatedFrame : public qglviewer::ManipulatedFrame
{
	public:
		void checkIfGrabsMouse(int x, int y, const qglviewer::Camera* const camera)
		{
			qglviewer::Vec proj = camera->projectedCoordinatesOf(position());
			setGrabsMouse((fabs(x-proj.x) < 50) && (fabs(y-proj.y) < 50)); // Rectangular region
		}

		void mousePressEvent(QMouseEvent *const event, qglviewer::Camera *const camera)
		{
			qglviewer::ManipulatedFrame::mousePressEvent(event, camera);

			moved=true;
		}

		void wheelEvent(QWheelEvent *const event, qglviewer::Camera *const camera)
		{
			qglviewer::ManipulatedFrame::wheelEvent(event, camera);

			moved=true;
		}

	public:
		bool moved;
};
typedef boost::shared_ptr<MeppManipulatedFrame> MeppManipulatedFramePtr;	//!< boost MeppManipulatedFrame shared pointeur

class mepp_component_plugin_interface;

/*! 
 * \class Viewer
 * \brief Viewer class.
 */
class Viewer : public QGLViewer
{
	Q_OBJECT

	public:
		/*!
		 * \brief Constructor.
		 *
		 * \param parent : parent window.
		 * \param lp : list of plugins/components.
		 */
		Viewer(QWidget *parent, QList<mepp_component_plugin_interface *> lp);
		/*!
		 * \brief Destructor.
		 */
		~Viewer();

		/*!
		 * \fn void WriteIni(bool force=false)
		 * \brief Write mepp.ini.
		 *
		 * \param force : force writing or not (default: false).
		 */
		void WriteIni(bool force=false);

		/*!
		 * \fn QString userFriendlyCurrentFile()
		 * \brief Return filename (without path) of the current active polyhedron.
		 *
		 * \return QString.
		 */
		QString userFriendlyCurrentFile() { return scene_ptr->userFriendlyCurrentFile(); }

		/*!
		 * \fn void setScenePtr(ScenePtr scene_ptr)
		 * \brief Set the Scene pointer to the viewer.
		 *
		 * \param scene_ptr the Scene pointer.
		 */
		void setScenePtr(ScenePtr scene_ptr);
		/*!
		 * \fn ScenePtr getScenePtr()
		 * \brief Get the Scene pointer of the viewer.
		 *
		 * \return the Scene pointer.
		 */
		ScenePtr getScenePtr() { return scene_ptr; }

		/*!
		 * \fn QWidget *getParent()
		 * \brief Get the parent widget of the viewer.
		 *
		 * \return QWidget*.
		 */
		QWidget *getParent() { return m_parent; }

		/*!
		 * \fn int get_nb_frames()
		 * \brief Get the number of frames (in Space mode).
		 *
		 * \return the number of frames.
		 */
		int get_nb_frames() { return (int)frame_.size(); }
		/*!
		 * \fn MeppManipulatedFrame* frame(unsigned short i)
		 * \brief Get the pointeur of the frame i (in Space mode).
		 *
		 * \param i the frame number.
		 * \return the pointeur of the frame.
		 */
		MeppManipulatedFrame* frame(unsigned short i) { return frame_[i].get(); }
		/*!
		 * \fn setSelectedFrameNumber(unsigned short nb)
		 * \brief Select a frame.
		 *
		 * \param nb the number of the frame to select.
		 */
		void setSelectedFrameNumber(unsigned short nb) { selected = nb; }
		/*!
		 * \fn addFrame()
		 * \brief Add an empty frame.
		 */
		void addFrame() { frame_.push_back(MeppManipulatedFramePtr(new MeppManipulatedFrame())); glList_.push_back(glGenLists(1)); }
		/*!
		 * \fn GLuint glList(unsigned short i)
		 * \brief Get the GL indice of the list i.
		 *
		 * \param i the number of the list.
		 * \return the GL indice of the list.
		 */
		GLuint glList(unsigned short i) { return glList_[i]; }

#if(0)
		/*!
		 * \fn qglviewer::Vec getInitialCameraPosition()
		 * \brief Get the camera position.
		 *
		 * \return the camera position.
		 */
		//qglviewer::Vec getInitialCameraPosition() { return initialCameraPosition; }
		/*!
		 * \fn qglviewer::Quaternion getInitialCameraOrientation()
		 * \brief Get the camera orientation.
		 *
		 * \return the camera orientation.
		 */
		//qglviewer::Quaternion getInitialCameraOrientation() { return initialCameraOrientation; }
#endif

		/*!
		 * \fn change_material(string mat_name)
		 * \brief Change the current material.
		 *
		 * \param mat_name the name of the current material (Silver, Gold, Jade, Light blue, Emerald, Polished silver, Chrome, Copper, Polished gold, Pewter, Obsidian, Black plastic, Polished bronze, Polished copper, Pearl, Ruby, Turquoise, Brass, None).
		 */
		//void change_material(string mat_name);
		#include "material.h"

		/*!
		 * \fn string get_material()
		 * \brief Get the current material.
		 *
		 * \return the name of the current material (Silver, Gold, Jade, Light blue, Emerald, Polished silver, Chrome, Copper, Polished gold, Pewter, Obsidian, Black plastic, Polished bronze, Polished copper, Pearl, Ruby, Turquoise, Brass, None).
		 */
		string get_material() { return m_last_material; }

		// rendering options
		/*!
		 * \fn setRender_Point()
		 * \brief Render mesh in point mode.
		 */
		void setRender_Point() { m_PolygonMode = GL_POINT; recreateListsAndUpdateGL(); }
		/*!
		 * \fn setRender_Line()
		 * \brief Render mesh in line mode.
		 */
		void setRender_Line() { m_PolygonMode = GL_LINE; recreateListsAndUpdateGL(); }
		/*!
		 * \fn setRender_Fill()
		 * \brief Render mesh in fill mode.
		 */
		void setRender_Fill() { m_PolygonMode = GL_FILL; recreateListsAndUpdateGL(); }
		/*!
		 * \fn int getRender_Mode()
		 * \brief Get the current render mode.
		 *
		 * \return the current render mode (GL_POINT, GL_LINE or GL_FILL).
		 */
		int getRender_Mode() { return m_PolygonMode; }

		/*!
		 * \fn setSuperimpose_Edges(bool b)
		 * \brief Superimpose (or not) edges of mesh(es).
		 *
		 * \param b true or false.
		 */
		void setSuperimpose_Edges(bool b) { m_SuperimposeEdges = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getSuperimpose_Edges
		 * \brief return true if superimpose edges of mesh(es) is active.
		 *
		 * \return true or false.
		 */
		bool getSuperimpose_Edges() { return m_SuperimposeEdges; }
		/*!
		 * \fn setSuperimpose_Vertices(bool b)
		 * \brief Superimpose (or not) vertices of mesh(es).
		 *
		 * \param b true or false.
		 */
		void setSuperimpose_Vertices(bool b) { m_SuperimposeVertices = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getSuperimpose_Vertices
		 * \brief return true if superimpose vertices of mesh(es) is active.
		 *
		 * \return true or false.
		 */
		bool getSuperimpose_Vertices() { return m_SuperimposeVertices; }
		/*!
		 * \fn setSuperimpose_Vertices_big(bool b)
		 * \brief Superimpose (or not) vertices (bigger mode) of mesh(es).
		 *
		 * \param b true or false.
		 */
		void setSuperimpose_Vertices_big(bool b) { m_SuperimposeVerticesBig = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getSuperimpose_Vertices_big
		 * \brief return true if superimpose vertices (bigger mode) of mesh(es) is active.
		 *
		 * \return true or false.
		 */
		bool getSuperimpose_Vertices_big() { return m_SuperimposeVerticesBig; }

		/*!
		 * \fn setVertex_Color(bool b)
		 * \brief View mesh(es) in vertex color mode (or not).
		 *
		 * \param b true or false.
		 */
		void setVertex_Color(bool b) { m_UseVertexColor = b; if (b) { m_UseFaceColor = !b; if (m_UseTexture) setTexture(false); } recreateListsAndUpdateGL(); }
		/*!
		 * \fn getVertex_Color
		 * \brief return true if viewing mesh(es) in vertex color mode is active.
		 *
		 * \return b true or false.
		 */
		bool getVertex_Color() { return m_UseVertexColor; }
		/*!
		 * \fn setFace_Color(bool b)
		 * \brief View mesh(es) in face color mode (or not).
		 *
		 * \param b true or false.
		 */
		void setFace_Color(bool b) { m_UseFaceColor = b; if (b) { m_UseVertexColor = !b; if (m_UseTexture) setTexture(false); } recreateListsAndUpdateGL(); }
		/*!
		 * \fn getFace_Color
		 * \brief return true if viewing mesh(es) in face color mode is active.
		 *
		 * \return b true or false.
		 */
		bool getFace_Color() { return m_UseFaceColor; }
		/*!
		 * \fn setTexture(bool b)
		 * \brief View mesh(es) in texture mode (or not).
		 *
		 * \param b true or false.
		 */
		void setTexture(bool b)
		{
			m_UseTexture = b;
			if (b)
			{
				m_last_material_saved = m_last_material;
				change_material("None");

				glEnable(GL_TEXTURE_2D);
				m_UseVertexColor = m_UseFaceColor = false;
			}
			else
			{
				m_last_material = m_last_material_saved;
				change_material(m_last_material);

				glDisable(GL_TEXTURE_2D);
			}
			recreateListsAndUpdateGL();
		}
		/*!
		 * \fn getTexture
		 * \brief return true if viewing mesh(es) in texture mode is active.
		 *
		 * \return b true or false.
		 */
		bool getTexture() { return m_UseTexture; }

		/*!
		 * \fn setLighting(bool b)
		 * \brief Toggle lighting (on/off).
		 *
		 * \param b true or false.
		 */
		void setLighting(bool b) { m_Lighting = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getLighting
		 * \brief return true if lighting is on.
		 *
		 * \return b true or false.
		 */
		bool getLighting() { return m_Lighting; }
		/*!
		 * \fn setSmooth_Shading(bool b)
		 * \brief Toggle smooth shading (on/off).
		 *
		 * \param b true or false.
		 */
		void setSmooth_Shading(bool b) { m_SmoothShading = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getSmooth_Shading
		 * \brief return true if smooth shading is on.
		 *
		 * \return b true or false.
		 */
		bool getSmooth_Shading() { return m_SmoothShading; }

		/*!
		 * \fn setAntialiasing(bool b)
		 * \brief Toggle antialiasing (on/off).
		 *
		 * \param b true or false.
		 */
		void setAntialiasing(bool b) { m_Antialiasing = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getAntialiasing
		 * \brief return true if antialiasing is on.
		 *
		 * \return b true or false.
		 */
		bool getAntialiasing() { return m_Antialiasing; }
		/*!
		 * \fn setCulling(bool b)
		 * \brief Toggle culling (on/off).
		 *
		 * \param b true or false.
		 */
		void setCulling(bool b) { m_Culling = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn getCulling
		 * \brief return true if culling is on.
		 *
		 * \return b true or false.
		 */
		bool getCulling() { return m_Culling; }
		// rendering options

		// color options
		/*!
		 * \fn setViewerBackgroundColor(QColor c)
		 * \brief Set background color.
		 *
		 * \param c a QColor.
		 */
		void setViewerBackgroundColor(QColor c) { m_BackColor[0] = float(c.red())/255.; m_BackColor[1] = float(c.green())/255.; m_BackColor[2] = float(c.blue())/255.; updateGL(); }
		/*!
		 * \fn QColor getViewerBackgroundColor()
		 * \brief Get background color.
		 *
		 * \return a QColor.
		 */
		QColor getViewerBackgroundColor() { return QColor(int(m_BackColor[0]*255.), int(m_BackColor[1]*255.), int(m_BackColor[2]*255.)); }

		/*!
		 * \fn setViewerVertexColor(QColor c)
		 * \brief Set vertex color.
		 *
		 * \param c a QColor.
		 */
		void setViewerVertexColor(QColor c) { m_VertexColor[0] = float(c.red())/255.; m_VertexColor[1] = float(c.green())/255.; m_VertexColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		/*!
		 * \fn QColor getViewerVertexColor()
		 * \brief Get vertex color.
		 *
		 * \return a QColor.
		 */
		QColor getViewerVertexColor() { return QColor(int(m_VertexColor[0]*255.), int(m_VertexColor[1]*255.), int(m_VertexColor[2]*255.)); }

		/*!
		 * \fn setViewerEdgeColor(QColor c)
		 * \brief Set edge color.
		 *
		 * \param c a QColor.
		 */
		void setViewerEdgeColor(QColor c) { m_EdgeColor[0] = float(c.red())/255.; m_EdgeColor[1] = float(c.green())/255.; m_EdgeColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		/*!
		 * \fn QColor getViewerEdgeColor()
		 * \brief Get edge color.
		 *
		 * \return a QColor.
		 */
		QColor getViewerEdgeColor() { return QColor(int(m_EdgeColor[0]*255.), int(m_EdgeColor[1]*255.), int(m_EdgeColor[2]*255.)); }

		/*!
		 * \fn setViewerFaceColor(QColor c)
		 * \brief Set face color.
		 *
		 * \param c a QColor.
		 */
		void setViewerFaceColor(QColor c) { m_MeshColor[0] = float(c.red())/255.; m_MeshColor[1] = float(c.green())/255.; m_MeshColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		/*!
		 * \fn QColor getViewerFaceColor()
		 * \brief Get face color.
		 *
		 * \return a QColor.
		 */
		QColor getViewerFaceColor() { return QColor(int(m_MeshColor[0]*255.), int(m_MeshColor[1]*255.), int(m_MeshColor[2]*255.)); }
		// color options

		// show options
		/*!
		 * \fn setShowNormals(bool b)
		 * \brief Toggle showing normals (on/off).
		 *
		 * \param b true or false.
		 */
		void setShowNormals(bool b) { show_normals = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn bool getShowNormals()
		 * \brief return true if normals are shown.
		 *
		 * \return b true or false.
		 */
		bool getShowNormals() { return show_normals; }

		/*!
		 * \fn setBounding_box(bool b)
		 * \brief Toggle showing bounding box (on/off).
		 *
		 * \param b true or false.
		 */
		void setBounding_box(bool b) { m_DrawBoundingBox = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn bool getBounding_box()
		 * \brief return true if bounding box is shown.
		 *
		 * \return b true or false.
		 */
		bool getBounding_box() { return m_DrawBoundingBox; }
		/*!
		 * \fn setBounding_box_when_moving(bool b)
		 * \brief Toggle showing bounding box when moving (on/off).
		 *
		 * \param b true or false.
		 */
		void setBounding_box_when_moving(bool b) { m_DrawBoundingBoxWhenMoving = b; recreateListsAndUpdateGL(); }
		/*!
		 * \fn bool getBounding_box_when_moving()
		 * \brief return true if bounding box when moving is shown.
		 *
		 * \return b true or false.
		 */
		bool getBounding_box_when_moving() { return m_DrawBoundingBoxWhenMoving; }
		// show options

		// view options
		/*!
		 * \fn showAllScene()
		 * \brief Set camera position and orientation to see the whole scene.
		 */
		void showAllScene()
		{
			setSceneRadius(1.0);
			setSceneCenter(qglviewer::Vec(0.0, 0.0, 0.0));
			camera()->setZNearCoefficient(0.005f);
			camera()->setZClippingCoefficient(sqrt(3.0));

			PolyhedronPtr p = scene_ptr->get_polyhedron();
			if (!p->empty())
			{
				setSceneBoundingBox(qglviewer::Vec(to_double(p->xmin()),to_double(p->ymin()),to_double(p->zmin())), qglviewer::Vec(to_double(p->xmax()),to_double(p->ymax()),to_double(p->zmax())));
				//camera()->setFOVToFitScene(); // MT
				showEntireScene();
			}

			p->pInitialCameraPosition = camera()->position();
			p->pInitialCameraOrientation = camera()->orientation();
		}

		/*!
		 * \fn showAllSceneForSpaceMode()
		 * \brief Set camera position and orientation to see the whole scene in Space mode.
		 */
		void showAllSceneForSpaceMode()
		{
			double xmin, ymin, zmin;
			double xmax, ymax, zmax;
			bool first = true;

			xmin = ymin = zmin = 0.;
			xmax = ymax = zmax = 0.;

			setSceneRadius(1.0);
			setSceneCenter(qglviewer::Vec(0.0, 0.0, 0.0));
			camera()->setZNearCoefficient(0.005f);
			camera()->setZClippingCoefficient(sqrt(3.0));

			int nbMesh = qMin(scene_ptr->get_nb_polyhedrons(), get_nb_frames());
			for (int i=0; i<nbMesh; i++)
			{
				if (!scene_ptr->get_polyhedron(i)->empty())
				{
					if (first)
					{
						xmin = to_double(scene_ptr->get_polyhedron(i)->xmin());
						ymin = to_double(scene_ptr->get_polyhedron(i)->ymin());
						zmin = to_double(scene_ptr->get_polyhedron(i)->zmin());

						xmax = to_double(scene_ptr->get_polyhedron(i)->xmax());
						ymax = to_double(scene_ptr->get_polyhedron(i)->ymax());
						zmax = to_double(scene_ptr->get_polyhedron(i)->zmax());

						first=false;
					}
					else
					{
						if (scene_ptr->get_polyhedron(i)->xmin() < xmin) xmin = scene_ptr->get_polyhedron(i)->xmin();
						if (scene_ptr->get_polyhedron(i)->ymin() < ymin) ymin = scene_ptr->get_polyhedron(i)->ymin();
						if (scene_ptr->get_polyhedron(i)->zmin() < zmin) zmin = scene_ptr->get_polyhedron(i)->zmin();

						if (scene_ptr->get_polyhedron(i)->xmax() > xmax) xmax = scene_ptr->get_polyhedron(i)->xmax();
						if (scene_ptr->get_polyhedron(i)->ymax() > ymax) ymax = scene_ptr->get_polyhedron(i)->ymax();
						if (scene_ptr->get_polyhedron(i)->zmax() > zmax) zmax = scene_ptr->get_polyhedron(i)->zmax();
					}
				}
			}
			
			if (!first)
			{
				setSceneBoundingBox(qglviewer::Vec(xmin,ymin,zmin), qglviewer::Vec(xmax,ymax,zmax));
				showEntireScene();
			}
		}

		/*!
		 * \fn centerAllObjects(bool forcePosition=true)
		 * \brief Center all objects in Space mode.
		 *
		 * \param forcePosition true or false (default: true).
		 */
		void centerAllObjects(bool forcePosition=true)
		{
			int nbMesh = qMin(scene_ptr->get_nb_polyhedrons(), get_nb_frames());
			for (int i=0; i<nbMesh; i++)
			{
				//if (!scene_ptr->get_polyhedron(i)->empty())
				{
					if (forcePosition)
						frame(i)->setPosition(qglviewer::Vec(0.f, 0.f, 0.f));
					else
					{
						int m = i%2;

						if (m)
							frame(i)->setPosition(qglviewer::Vec(0.f, getYStep() * ceil(i/2.), 0.f));
						else
							frame(i)->setPosition(qglviewer::Vec(0.f, -getYStep() * ceil(i/2.), 0.f));
					}

					frame(i)->setOrientation(qglviewer::Quaternion(0, 0, 0, 1)); // identity Quaternion

					frame(i)->moved=true;
				}
			}

			updateGL();
		}

#if(0)
		/*!
		 * \fn resetView()
		 * \brief Reset the camera position and orientation to see the initial view.
		 */
		void resetView()
		{
			/*camera()->setPosition(getInitialCameraPosition());
			camera()->setOrientation(getInitialCameraOrientation());*/

			PolyhedronPtr p = scene_ptr->get_polyhedron();
			//if (!p->empty())
			{
				camera()->setPosition(p->pInitialCameraPosition);
				camera()->setOrientation(p->pInitialCameraOrientation);
			}

			updateGL();
		}
#endif

		/*!
		 * \fn setCouplingTranslations(bool b)
		 * \brief Coupling translations of meshes (on/off).
		 *
		 * \param b true or false.
		 */
		void setCouplingTranslations(bool b) { mCouplingTranslations = b; updateGL(); }
		/*!
		 * \fn bool getCouplingTranslations()
		 * \brief return true if coupling translations of meshes is on.
		 *
		 * \return b true or false.
		 */
		bool getCouplingTranslations() { return mCouplingTranslations; }

		/*!
		 * \fn setCouplingRotations(bool b)
		 * \brief Coupling rotations of meshes (on/off).
		 *
		 * \param b true or false.
		 */
		void setCouplingRotations(bool b) { mCouplingRotations = b; updateGL(); }
		/*!
		 * \fn bool getCouplingRotations()
		 * \brief return true if coupling rotations of meshes is on.
		 *
		 * \return b true or false.
		 */
		bool getCouplingRotations() { return mCouplingRotations; }

		/*!
		 * \fn setCouplingZooms(bool b)
		 * \brief Coupling zooms of meshes (on/off).
		 *
		 * \param b true or false.
		 */
		void setCouplingZooms(bool b) { mCouplingZooms = b; updateGL(); }
		/*!
		 * \fn bool getCouplingZooms()
		 * \brief return true if coupling zooms of meshes is on.
		 *
		 * \return b true or false.
		 */
		bool getCouplingZooms() { return mCouplingZooms; }

		/*!
		 * \fn setVBO_mode(bool b)
		 * \brief Toggle display lists mode (on/off).
		 *
		 * \param b true or false.
		 */
        void setVBO_mode(bool b) { VBO_mode = b; if (b) setVBO_modeUncheck(b); /*setMouseTracking(!b);*/ recreateListsAndUpdateGL(); } // pb with setMouseTracking
		/*!
		 * \fn bool getVBO_mode()
		 * \brief return true if display lists mode is on.
		 *
		 * \return b true or false.
		 */
		bool getVBO_mode() { return VBO_mode; }
		// view options

		// capture options
		/*!
		 * \fn setSave_animation(bool b)
		 * \brief Toggle the fact that we are currently saving an animation or not (on/off).
		 *
		 * \param b true or false.
		 */
		void setSave_animation(bool b) { save_animation = b; }
		/*!
		 * \fn bool getSave_animation()
		 * \brief return true if we are currently saving an animation.
		 *
		 * \return b true or false.
		 */
		bool getSave_animation() { return save_animation; }

		/*!
		 * \fn saveAnimation(bool b, QString &savePNGLocation)
		 * \brief Start or stop saving an animation (on/off).
		 *
		 * \param b true (on) or false (off).
		 * \param savePNGLocation contains the default PNG save location.
		 */
		void saveAnimation(bool b, QString &savePNGLocation)
		{
			if (b)
			{
				QString fileName = QFileDialog::getSaveFileName(this, QObject::tr("Save Screenshot Sequence"),
										 /*QDir::currentPath()*/savePNGLocation,
										 QObject::tr("PNG Files (*.png)"));
				if (!fileName.isEmpty())
				{
					savePNGLocation = QFileInfo(fileName).absolutePath();

					QFileInfo fileInfo(fileName);
					fileName = fileInfo.absolutePath()+ '/' + fileInfo.baseName() + ".png";

					setSnapshotFileName(fileName);
					setSnapshotFormat("PNG");

					saveSnapshot(true, true);

					connect(this, SIGNAL(drawFinished(bool)), SLOT(saveSnapshot(bool)));
					setSave_animation(true);
				}
				else
					setSave_animation(false);
			}
			else
			{
				disconnect(SIGNAL(drawFinished(bool)));
				setSave_animation(false);
			}
		}

#ifdef WITH_FFMPEG
		/*!
		 * \fn saveFFmpegAnimation(bool b, QString &saveAVILocation, int bitrate=1500, int fps=25, bool resize=false, int gop=12)
		 * \brief Start or stop saving an FFmpeg animation (on/off).
		 *
		 * \param b true or false.
		 * \param saveAVILocation contains the default AVI save location.
		 * \param bitrate bitrate of the video (default: 1500).
                 * \param fps fps of the video (default: 25).
 		 * \param resize false (unused).
 		 * \param gop gop of the video (default: 12).                 
		 */
		void saveFFmpegAnimation(bool b, QString &saveAVILocation, int bitrate=1500, int fps=25, bool resize=false, int gop=12)
		{
			if (b)
			{
				QString fileName = QFileDialog::getSaveFileName(this, QObject::tr("Save MPEG-4 Video"),
										 /*QDir::currentPath()*/saveAVILocation,
										 QObject::tr("AVI Files (*.avi)"));
				if (!fileName.isEmpty())
				{
					saveAVILocation = QFileInfo(fileName).absolutePath();

					QFileInfo fileInfo(fileName);
					fileName = fileInfo.absolutePath()+ '/' + fileInfo.baseName() + ".avi";

					/*if (resize)
						this->resize(720, 576); // PAL DV format*/

					int width=this->width();
					int height=this->height();
					//int bitrate=1000;
					//int gop=20;
					encoder.createFile(fileName,width,height,bitrate*1000,gop,fps);

					connect(this, SIGNAL(drawFinished(bool)), SLOT(shotCapture(bool)));
					setSave_animation(true);
				}
				else
					setSave_animation(false);
			}
			else
			{
				disconnect(SIGNAL(drawFinished(bool)));
				setSave_animation(false);

				encoder.close();
			}
		}
#endif
		// capture options

		// dynamic options
		/*!
		 * \fn setDynTitle()
		 * \brief set and update the title of the viewer window.
		 */
		void setDynTitle()
		{
			if (getScenePtr()->get_loadType() == Normal)
				setWindowTitle(QObject::tr("vid: %2 - %1")
										.arg(userFriendlyCurrentFile())
										.arg((qlonglong)this, 0, 16));
			else
				setWindowTitle(QObject::tr("vid: %2 - (%3: %4/%5) - %1")
										.arg(userFriendlyCurrentFile())
										.arg((qlonglong)this, 0, 16)
										.arg(getScenePtr()->get_stringLoadType())
										.arg(getScenePtr()->get_current_polyhedron()+1)
										.arg(getScenePtr()->get_nb_polyhedrons()));
		}

		/*!
		 * \fn setFps(int fps)
		 * \brief set fps before playing dynamic time sequence.
		 *
		 * \param fps the fps.
		 */
		void setFps(int fps) { m_fps = fps; }
		/*!
		 * \fn int getFps()
		 * \brief get the fps used by dynamic time sequence.
		 *
		 * \return the fps.
		 */
		int getFps() { return m_fps; }

		/*!
		 * \fn setDynReverseStart()
		 * \brief Start playing dynamic time sequence in reverse order.
		 */
		void setDynReverseStart() { m_reverse = true; m_loop = false; timerDynamic->start(1000/m_fps); }
		/*!
		 * \fn setDynReverseStartLoop()
		 * \brief Start playing dynamic time sequence in reverse order with loop.
		 */
		void setDynReverseStartLoop() { m_reverse = true; m_loop = true; timerDynamic->start(1000/m_fps); }
		/*!
		 * \fn setDynStart()
		 * \brief Start playing dynamic time sequence.
		 */
		void setDynStart() { m_reverse = false; m_loop = false; timerDynamic->start(1000/m_fps); }
		/*!
		 * \fn setDynStartLoop()
		 * \brief Start playing dynamic time sequence with loop.
		 */
		void setDynStartLoop() { m_reverse = false; m_loop = true; timerDynamic->start(1000/m_fps); }
		/*!
		 * \fn setDynStop()
		 * \brief Stop playing dynamic time sequence.
		 */
		void setDynStop() { timerDynamic->stop(); }

		/*!
		 * \fn setDynFirst()
		 * \brief Go to first position of dynamic time sequence.
		 */
		void setDynFirst() { scene_ptr->set_current_polyhedron(0); setDynTitle(); recreateListsAndUpdateGL(); }
		/*!
		 * \fn setDynPrevious()
		 * \brief Go to previous position of dynamic time sequence.
		 */
		void setDynPrevious() { if (scene_ptr->get_current_polyhedron() >= 1) scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()-1); setDynTitle(); recreateListsAndUpdateGL(); }
		/*!
		 * \fn setDynNext()
		 * \brief Go to next position of dynamic time sequence.
		 */
		void setDynNext() { if (scene_ptr->get_current_polyhedron() < (scene_ptr->get_nb_polyhedrons()-1)) scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()+1); setDynTitle(); recreateListsAndUpdateGL(); }
		/*!
		 * \fn setDynLast()
		 * \brief Go to last position of dynamic time sequence.
		 */
		void setDynLast() { scene_ptr->set_current_polyhedron(scene_ptr->get_nb_polyhedrons()-1); setDynTitle(); recreateListsAndUpdateGL(); }
		// dynamic options

		/*!
		 * \fn setDelete()
		 * \brief Delete current mesh - only in Space or Time mode.
		 */
		void setDelete()
		{
			if (scene_ptr->get_nb_polyhedrons() > 1)
			{
				int i = scene_ptr->get_current_polyhedron();
				int id = i;

				scene_ptr->delete_polyhedron(i);
				if (i!=0)
					i=i-1;
				scene_ptr->set_current_polyhedron(i);

				// for Space mode consistency
				for (int f=id; f<get_nb_frames(); f++)
				{
					if ((f+1)<get_nb_frames())
					{
						frame(f)->setPosition(frame(f+1)->position());
						frame(f)->setOrientation(frame(f+1)->orientation());
					}
				}
				for (int f=0; f<get_nb_frames(); f++)
				{
					if (f > scene_ptr->get_nb_polyhedrons()-1)
					{
						int m = f%2;
						if (m)
							frame(f)->setPosition(qglviewer::Vec(0.f, mYStep * ceil(f/2.), 0.f));
						else
							frame(f)->setPosition(qglviewer::Vec(0.f, -mYStep * ceil(f/2.), 0.f));

						frame(f)->setOrientation(qglviewer::Quaternion(0, 0, 0, 1)); // identity Quaternion
					}
				}
				// for Space mode consistency

				setDynTitle();
				recreateListsAndUpdateGL();
			}
			else
				QMessageBox::warning(m_parent, APPLICATION, QObject::tr("Deleting last mesh not allowed."));
		}

		/*!
		 * \fn recreateListsAndUpdateGL()
		 * \brief Recreate all GL lists (if display lists are active) and update GL.
		 */
		void recreateListsAndUpdateGL()
		{
			if (VBO_mode)
				createLists = true;

			updateGL();

			/*if (VBO_mode)
				createLists = false;*/
		}

		/*!
		 * \fn double getYStep()
		 * \brief get the 'y step' used for Space Mode (between each mesh).
		 *
		 * \return the 'y step'.
		 */
		double getYStep() { return mYStep; }

		/*!
		 * \fn bool getForceFFMPEG()
		 * \brief return true if we want always FFMPEG if available.
		 *
		 * \return true of false.
		 */
		bool getForceFFMPEG() { return mForceFFMPEGifAvailable; }

	protected:
		/*!
		 * \fn init()
		 * \brief Init of the viewer.
		 */
		virtual void init();
		/*!
		 * \fn postSelection(const QPoint& point)
		 * \brief Function to select a specific frame with mouse (Space mode).
		 *
		 * \param point the point in screen coordinates.
		 */
		virtual void postSelection(const QPoint& point);
		/*!
		 * \fn drawWithNames()
		 * \brief Drawing (with names) of the viewer, used for selecting in Space mode.
		 */
		virtual void drawWithNames();
		/*!
		 * \fn draw()
		 * \brief Drawing of the viewer.
		 */
		virtual void draw();
		/*!
		 * \fn helpString()
		 * \brief Function to show an help for QGLViewer.
		 */
		virtual QString helpString() const;

		/*!
		 * \fn render(bool sel, bool grab)
		 * \brief Main render function.
		 *
		 * \param sel true if the current frame is selected. Always false if Normal or Time mode.
		 * \param grab true if the mouse is over the current frame. Always false if Normal or Time mode.
		 */
		void render(bool sel, bool grab);

		/*!
		 * \brief Close the viewer (event).
		 *
		 * \param event the event.
		 */
		void closeEvent(QCloseEvent *event);

		/*!
		 * \brief Show the context menu (event). Not used.
		 *
		 * \param event the event.
		 */
		void contextMenuEvent(QContextMenuEvent *event);
		/*!
		 * \brief Show the context menu (event).
		 *
		 * \param event the event.
		 */
		void MEPPcontextMenuEvent(QMouseEvent *event);

		/*!
		 * \fn setVBO_modeUncheck(bool b)
		 * \brief Used to disable or not the bounding box moving if display lists are on or not.
		 *
		 * \param b true or false.
		 */
		void setVBO_modeUncheck(bool b);

		// events
		/*!
		 * \brief Event when a bouton button is pressed (event).
		 *
		 * \param event the event.
		 */
		virtual void mousePressEvent(QMouseEvent *event);
		/*!
		 * \brief Event when mouse is moved (event).
		 *
		 * \param event the event.
		 */
		virtual void mouseMoveEvent(QMouseEvent *event);
		/*!
		 * \brief Event when a bouton button is released (event).
		 *
		 * \param event the event.
		 */
		virtual void mouseReleaseEvent(QMouseEvent *event);
		/*!
		 * \brief Event when mouse wheel is moved (event).
		 *
		 * \param event the event.
		 */
		virtual void wheelEvent(QWheelEvent *event);
		/*!
		 * \brief Event when a key is pressed (event).
		 *
		 * \param event the event.
		 */
		virtual void keyPressEvent(QKeyEvent *event);
		/*!
		 * \brief Event when a key is released (event).
		 *
		 * \param event the event.
		 */
		virtual void keyReleaseEvent(QKeyEvent *event);

	private slots:
		/*!
		 * \fn shotDynamic
		 * \brief Used to update dynamic time sequence.
		 */
		void shotDynamic()
		{
			if (!m_reverse)
			{
				if (scene_ptr->get_current_polyhedron() < (scene_ptr->get_nb_polyhedrons()-1))
					scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()+1);
				else
				{
					if (m_loop)
						scene_ptr->set_current_polyhedron(0);
					else
					{
						timerDynamic->stop();
						return;
					}
				}
			}
			else
			{
				if (scene_ptr->get_current_polyhedron() > 0)
					scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()-1);
				else
				{
					if (m_loop)
						scene_ptr->set_current_polyhedron((scene_ptr->get_nb_polyhedrons()-1));
					else
					{
						timerDynamic->stop();
						return;
					}
				}
			}

			setDynTitle();
			recreateListsAndUpdateGL();
		}


		/*!
		 * \fn shotCapture
		 * \brief Used to encode FFmpeg animation.
		 *
		 * \param param not used, always true.
		 */
		void shotCapture(bool param = true)
		{
#ifdef WITH_FFMPEG
			this->makeCurrent();
			this->raise();
			encoder.encodeImage(this->grabFrameBuffer(true));
#endif
		}

		/*!
		 * \fn setActivePolyhedron(int p)
		 * \brief Set the current active polyhedron.
		 *
		 * \param p the polyhedron.
		 */
		void setActivePolyhedron(int p);

	public:
		QList<mepp_component_plugin_interface *> lplugin;	//!< list of plugins/components

	private:
		QWidget *m_parent;			//!< parent widget of the viewer
		ScenePtr scene_ptr;			//!< Scene pointer associated to the viewer

		bool m_AutoSaveIni;

		// OpenGL
		bool m_Lighting;
		bool m_Culling;
		int m_PolygonMode;
		bool m_SmoothShading;
		bool m_UseVertexColor;
		bool m_UseFaceColor;
		bool m_UseTexture;
		bool m_UseNormals;
		//bool m_FirstView;
		bool m_SuperimposeEdges;
		bool m_SuperimposeVertices;
		bool m_SuperimposeVerticesBig;
		bool m_Antialiasing;
		//float m_ThicknessControlEdges;
		//float m_PointSize;
		bool m_DrawBoundingBox;
		bool m_DrawBoundingBoxWhenMoving;
		bool m_DrawVoronoiEdges;

		// colors
		float m_BackColor[3];
		float m_MeshColor[3];
		float m_EdgeColor[3];
		float m_VertexColor[3];

		// mouse
		bool m_LeftButtonDown;
		bool m_RightButtonDown;
		bool m_Moving;
		bool m_HasMoved;

		string m_last_material, m_last_material_saved;

		//

		/*qglviewer::Vec initialCameraPosition;
		qglviewer::Quaternion initialCameraOrientation;*/

		//

		vector<MeppManipulatedFramePtr> frame_;
		vector<GLuint> glList_;
		unsigned short selected;

		/*!
		 * \fn dessine_space(bool names=false)
		 * \brief Used to draw the scene in Space mode.
		 *
		 * \param names if true we are drawing with names which is used for selecting in Space mode.
		 */
		void dessine_space(bool names=false);

		/*!
		 * \fn dessine(bool names=false)
		 * \brief Used to draw the scene in Normal and Time mode.
		 *
		 * \param names always false.
		 */
		void dessine(bool names=false);

		//

		bool mCouplingTranslations;
		bool mCouplingRotations;
		bool mCouplingZooms;

		double mYStep;

		bool show_normals;
		bool VBO_mode;
		bool save_animation;
		bool mForceFFMPEGifAvailable;

		//

		QTimer *timerDynamic;
		int m_fps;
		bool m_reverse, m_loop;

		//

		GLuint glId;
		bool createLists;

		//

		qglviewer::Vec orig, dir, selectedPoint;

#ifdef WITH_FFMPEG
		QVideoEncoder encoder;
#endif
};

#endif
