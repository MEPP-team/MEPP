///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 

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
		if (event->mouseButtons() & Qt::RightButton) // because Mac OS X
			bType=bRight;
		else if (event->mouseButtons() & Qt::LeftButton)
			bType=bLeft; 
		else
			bType=bNone;

		event->acceptProposedAction();
	}
}

void MdiArea::dropEvent(QDropEvent *event)
{
	int res = 0;
	Viewer *viewer = NULL;
	QList<QUrl> urls = event->mimeData()->urls();

	for (int i=0; i<urls.size(); i++)
	{
		QFileInfo fi(urls[i].toLocalFile());
		QString ext = fi.suffix();

		if ((ext.toLower()=="off") || (ext.toLower()=="obj") || (ext.toLower()=="smf") || (ext.toLower()=="ply"))
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