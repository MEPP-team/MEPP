/*!
 * \file mdiarea.cpp
 * \brief MdiArea file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */
#include "mdiarea.hxx"
#include "mainwindow.hxx"

MdiArea::MdiArea(QWidget *parent) : QMdiArea(parent)
{
	bType = bNone;
}

void MdiArea::paintEvent(QPaintEvent *paintEvent)
{
	QMdiArea::paintEvent(paintEvent);

	// and now only paint your image here
	QPainter painter(viewport());
 
	painter.fillRect(paintEvent->rect(), QColor(23, 74, 124));
	//painter.drawImage(paintEvent->rect()/*QPoint(0, 0)*/, QImage("./mepp_background.bmp"));
 
	painter.end();
}

void MdiArea::dragEnterEvent(QDragEnterEvent *event)
{
	if (event->mimeData()->hasUrls())	//hasFormat("text/uri-list")
	{
		event->acceptProposedAction();

// ---------------------------------------------------
#ifdef __linux__
	bType=bLeft;
	if (event->mouseButtons() & Qt::RightButton)
		bType=bRight;
#endif
// ---------------------------------------------------
	}
}

void MdiArea::dropEvent(QDropEvent *event)
{
	int res = 0;
	Viewer *viewer = NULL;
	QList<QUrl> urls = event->mimeData()->urls();

// ---------------------------------------------------
#ifdef __APPLE__
	bType=bLeft;
	if (event->keyboardModifiers() & Qt::MetaModifier)
		bType=bRight;
#endif
#ifdef _MSC_VER
	bType=bLeft;
	if (event->mouseButtons() & Qt::RightButton)
		bType=bRight;
#endif
// ---------------------------------------------------

	for (int i=0; i<urls.size(); i++)
	{
		QFileInfo fi(urls[i].toLocalFile());
		QString ext = fi.suffix();

		if ((ext.toLower()=="off") || (ext.toLower()=="obj") || (ext.toLower()=="smf") || (ext.toLower()=="ply") || (ext.toLower()=="x3d"))
		{
			if (m_mw->activeMdiChild() == 0)
				res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
			else
			{
				viewer = qobject_cast<Viewer *>(m_mw->activeMdiChild());

				if (bType == bLeft)
					res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
				else if (bType == bRight)
				{
	#ifdef __linux__
					if (event->keyboardModifiers() & Qt::MetaModifier)
	#else
					if (event->keyboardModifiers() & Qt::AltModifier)
	#endif
						res = m_mw->addFile(viewer, urls[i].toLocalFile(), Time, NULL);
					else
						res = m_mw->addFile(viewer, urls[i].toLocalFile(), Space, NULL);
				}
			}

			if (res)
				break;
		}
	}
	if (viewer)
		viewer->recreateListsAndUpdateGL();

	event->acceptProposedAction();
}