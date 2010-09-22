///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_VIEWER
#define HEADER_VIEWER

#include <QGLViewer/qglviewer.h>

#include <QtGui/QCloseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>

#include "scene.h"

typedef boost::shared_ptr<qglviewer::ManipulatedFrame> ManipulatedFramePtr;

class mepp_component_plugin_interface;

class Viewer : public QGLViewer
{
	Q_OBJECT

	public:
		Viewer(QWidget *parent, QList<mepp_component_plugin_interface *> lp);
		~Viewer();

		QString userFriendlyCurrentFile() { return scene_ptr->userFriendlyCurrentFile(); }

		void setScenePtr(ScenePtr scene_ptr);
		ScenePtr getScenePtr() { return scene_ptr; }

		QWidget *getParent() { return m_parent; }

		int get_nb_frames() { return (int)frame_.size(); }
		qglviewer::ManipulatedFrame* frame(unsigned short i) { return frame_[i].get(); }
		void setSelectedFrameNumber(unsigned short nb) { selected = nb; }
		void addFrame() { frame_.push_back(ManipulatedFramePtr(new qglviewer::ManipulatedFrame())); glList_.push_back(glGenLists(1)); if (VBO_mode) createLists = true; }
		GLuint glList(unsigned short i) { return glList_[i]; }

		qglviewer::Vec getInitialCameraPosition() { return initialCameraPosition; }
		qglviewer::Quaternion getInitialCameraOrientation() { return initialCameraOrientation; }

		void change_material(string mat_name);
		string get_material() { return m_last_material; }

		// rendering options
		void setRender_Point() { m_PolygonMode = GL_POINT; recreateListsAndUpdateGL(); }
		void setRender_Line() { m_PolygonMode = GL_LINE; recreateListsAndUpdateGL(); }
		void setRender_Fill() { m_PolygonMode = GL_FILL; recreateListsAndUpdateGL(); }	
		int getRender_Mode() { return m_PolygonMode; }

		void setSuperimpose_Edges(bool b) { m_SuperimposeEdges = b; recreateListsAndUpdateGL(); }
		bool getSuperimpose_Edges() { return m_SuperimposeEdges; }
		void setSuperimpose_Vertices(bool b) { m_SuperimposeVertices = b; recreateListsAndUpdateGL(); }
		bool getSuperimpose_Vertices() { return m_SuperimposeVertices; }
		void setSuperimpose_Vertices_big(bool b) { m_SuperimposeVerticesBig = b; recreateListsAndUpdateGL(); }
		bool getSuperimpose_Vertices_big() { return m_SuperimposeVerticesBig; }

		void setVertex_Color(bool b) { m_UseVertexColor = b; if (b) m_UseFaceColor = !b; recreateListsAndUpdateGL(); }
		bool getVertex_Color() { return m_UseVertexColor; }
		void setFace_Color(bool b) { m_UseFaceColor = b; if (b) m_UseVertexColor = !b; recreateListsAndUpdateGL(); }
		bool getFace_Color() { return m_UseFaceColor; }

		void setLighting(bool b) { m_Lighting = b; recreateListsAndUpdateGL(); }
		bool getLighting() { return m_Lighting; }
		void setSmooth_Shading(bool b) { m_SmoothShading = b; recreateListsAndUpdateGL(); }
		bool getSmooth_Shading() { return m_SmoothShading; }

		void setAntialiasing(bool b) { m_Antialiasing = b; recreateListsAndUpdateGL(); }
		bool getAntialiasing() { return m_Antialiasing; }
		void setCulling(bool b) { m_Culling = b; recreateListsAndUpdateGL(); }
		bool getCulling() { return m_Culling; }
		// rendering options

		// color options
		void setViewerBackgroundColor(QColor c) { m_BackColor[0] = float(c.red())/255.; m_BackColor[1] = float(c.green())/255.; m_BackColor[2] = float(c.blue())/255.; updateGL(); }
		QColor getViewerBackgroundColor() { return QColor(int(m_BackColor[0]*255.), int(m_BackColor[1]*255.), int(m_BackColor[2]*255.)); }

		void setViewerVertexColor(QColor c) { m_VertexColor[0] = float(c.red())/255.; m_VertexColor[1] = float(c.green())/255.; m_VertexColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		QColor getViewerVertexColor() { return QColor(int(m_VertexColor[0]*255.), int(m_VertexColor[1]*255.), int(m_VertexColor[2]*255.)); }

		void setViewerEdgeColor(QColor c) { m_EdgeColor[0] = float(c.red())/255.; m_EdgeColor[1] = float(c.green())/255.; m_EdgeColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		QColor getViewerEdgeColor() { return QColor(int(m_EdgeColor[0]*255.), int(m_EdgeColor[1]*255.), int(m_EdgeColor[2]*255.)); }

		void setViewerFaceColor(QColor c) { m_MeshColor[0] = float(c.red())/255.; m_MeshColor[1] = float(c.green())/255.; m_MeshColor[2] = float(c.blue())/255.; recreateListsAndUpdateGL(); }
		QColor getViewerFaceColor() { return QColor(int(m_MeshColor[0]*255.), int(m_MeshColor[1]*255.), int(m_MeshColor[2]*255.)); }
		// color options

		// show options
		void setShowNormals(bool b) { show_normals = b; recreateListsAndUpdateGL(); }
		bool getShowNormals() { return show_normals; }

		void setBounding_box(bool b) { m_DrawBoundingBox = b; recreateListsAndUpdateGL(); }
		bool getBounding_box() { return m_DrawBoundingBox; }
		void setBounding_box_when_moving(bool b) { m_DrawBoundingBoxWhenMoving = b; recreateListsAndUpdateGL(); }
		bool getBounding_box_when_moving() { return m_DrawBoundingBoxWhenMoving; }
		// show options

		// view options
		void showAllScene()
		{
			PolyhedronPtr p = scene_ptr->get_polyhedron();

			if (!p->empty())
			{
				setSceneBoundingBox(qglviewer::Vec(p->xmin(),p->ymin(),p->zmin()), qglviewer::Vec(p->xmax(),p->ymax(),p->zmax()));		
				showEntireScene();
			}
		}

		void centerAllObjects()
		{
			int nbMesh = qMin(scene_ptr->get_nb_polyhedrons(), get_nb_frames());
			for (int i=0; i<nbMesh; i++)
			{
				if (!scene_ptr->get_polyhedron(i)->empty())
				{
					frame(i)->setPosition(qglviewer::Vec(0.f, 0.f, 0.f));
					frame(i)->setOrientation(qglviewer::Quaternion(0, 0, 0, 1)); // identity Quaternion
				}
			}

			updateGL();
		}

		void resetView()
		{
			camera()->setPosition(getInitialCameraPosition());
			camera()->setOrientation(getInitialCameraOrientation());

			updateGL();
		}

		void setCouplingTranslations(bool b) { mCouplingTranslations = b; updateGL(); }
		bool getCouplingTranslations(bool b) { return mCouplingTranslations; }

		void setCouplingRotations(bool b) { mCouplingRotations = b; updateGL(); }
		bool getCouplingRotations() { return mCouplingRotations; }

		void setCouplingZooms(bool b) { mCouplingZooms = b; updateGL(); }
		bool getCouplingZooms(bool b) { return mCouplingZooms; }

                void setVBO_mode(bool b) { VBO_mode = b; if (b) setVBO_modeUncheck(b); recreateListsAndUpdateGL(); }
		bool getVBO_mode() { return VBO_mode; }
		// view options

		// capture options
		void setSave_animation(bool b) { save_animation = b; }
		bool getSave_animation() { return save_animation; }

		void saveAnimation(bool b)
		{
			if (b)
			{
				QString fileName = QFileDialog::getSaveFileName(this, tr("Save Screenshot Sequence"),
										 QDir::currentPath(),
										 tr("PNG Files (*.png)"));
				if (!fileName.isEmpty())
				{
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
		// capture options

		// dynamic options
		void setDynTitle() { setWindowTitle(tr("%1 - (%2: %3/%4)")
										.arg(userFriendlyCurrentFile())
										.arg(getScenePtr()->get_stringLoadType())
										.arg(getScenePtr()->get_current_polyhedron()+1)
										.arg(getScenePtr()->get_nb_polyhedrons())); }

		void setDynParams() { QMessageBox::information(m_parent, APPLICATION, tr("Function not yet implemented.")); }

		void setDynReverseStart() { m_reverse = true; m_loop = false; timerDynamic->start(1000/m_fps); }
		void setDynReverseStartLoop() { m_reverse = true; m_loop = true; timerDynamic->start(1000/m_fps); }
		void setDynStart() { m_reverse = false; m_loop = false; timerDynamic->start(1000/m_fps); }
		void setDynStartLoop() { m_reverse = false; m_loop = true; timerDynamic->start(1000/m_fps); }
		void setDynStop() { timerDynamic->stop(); }

		void setDynFirst() { scene_ptr->set_current_polyhedron(0); setDynTitle(); recreateListsAndUpdateGL(); }
		void setDynPrevious() { if (scene_ptr->get_current_polyhedron() >= 1) scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()-1); setDynTitle(); recreateListsAndUpdateGL(); }
		void setDynNext() { if (scene_ptr->get_current_polyhedron() < (scene_ptr->get_nb_polyhedrons()-1)) scene_ptr->set_current_polyhedron(scene_ptr->get_current_polyhedron()+1); setDynTitle(); recreateListsAndUpdateGL(); }
		void setDynLast() { scene_ptr->set_current_polyhedron(scene_ptr->get_nb_polyhedrons()-1); setDynTitle(); recreateListsAndUpdateGL(); }

		void setDynDelete()
		{
			if (scene_ptr->get_nb_polyhedrons() > 1)
			{
				int i = scene_ptr->get_current_polyhedron();

				scene_ptr->delete_polyhedron(i);
				if (i!=0)
					i=i--;
				scene_ptr->set_current_polyhedron(i);

				setDynTitle();
				recreateListsAndUpdateGL();
			}
			else
				QMessageBox::warning(m_parent, APPLICATION, tr("Deleting last mesh not allowed."));
		}
		// dynamic options

		void recreateListsAndUpdateGL()
		{
			if (VBO_mode)
				createLists = true;

			updateGL();
		}

	protected:
		virtual void init();
		virtual void postSelection(const QPoint& point);
		virtual void drawWithNames();
		virtual void draw();
		virtual QString helpString() const;

		void render();

		void closeEvent(QCloseEvent *event);

                void setVBO_modeUncheck(bool b);

		// events
		virtual void mousePressEvent(QMouseEvent *event);
		virtual void mouseMoveEvent(QMouseEvent *event);
		virtual void mouseReleaseEvent(QMouseEvent *event);
		virtual void wheelEvent(QWheelEvent *event);
		virtual void keyPressEvent(QKeyEvent *event);
		virtual void keyReleaseEvent(QKeyEvent *event);

	private slots:
		void shotDynamic();
		void shotCapture();

	public:
		QList<mepp_component_plugin_interface *> lplugin;

	private:
		QWidget *m_parent;
		ScenePtr scene_ptr;

		// OpenGL
		bool m_Lighting;
		bool m_Culling;
		int m_PolygonMode;
		bool m_SmoothShading;
		bool m_UseVertexColor;
		bool m_UseFaceColor;
		bool m_UseNormals;
		bool m_FirstView;
		bool m_SuperimposeEdges;
		bool m_SuperimposeVertices;
		bool m_SuperimposeVerticesBig;
		bool m_Antialiasing;
		float m_ThicknessControlEdges;
		float m_PointSize;
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

		string m_last_material;

		//

		qglviewer::Vec initialCameraPosition;
		qglviewer::Quaternion initialCameraOrientation;

		//

		vector<ManipulatedFramePtr> frame_;
		vector<GLuint> glList_;
		unsigned short selected;

		void dessine_space(bool names=false);

		void dessine(bool names=false);

		//

		bool mCouplingTranslations;
		bool mCouplingRotations;
		bool mCouplingZooms;

		bool show_normals;
		bool VBO_mode;
		bool save_animation;

		//

		QTimer *timerDynamic;
		int m_fps;
		bool m_reverse, m_loop;

		//

		GLuint glId;
		bool createLists;
};

#endif
