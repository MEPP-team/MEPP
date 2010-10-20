///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 

#include "mdiarea.hxx"
#include "mainwindow.hxx"

MdiArea::MdiArea(QWidget *parent) : QMdiArea(parent)
{
}

void MdiArea::dragEnterEvent(QDragEnterEvent *event)
{
	if (event->mimeData()->hasUrls())	//hasFormat("text/uri-list")
		event->acceptProposedAction();
}

void MdiArea::dropEvent(QDropEvent *event)
{
	int res = 0;
	Viewer *viewer = NULL;
	QList<QUrl> urls = event->mimeData()->urls();

	for (int i=0; i<urls.size(); i++)
	{
		QFileInfo fi(urls[i].toLocalFile());
		QString ext = fi.completeSuffix();

		if ((ext.toLower()=="off") || (ext.toLower()=="obj") || (ext.toLower()=="smf") || (ext.toLower()=="ply"))
		{
			if (m_mw->activeMdiChild() == 0)
				res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
			else
			{
				viewer = qobject_cast<Viewer *>(m_mw->activeMdiChild());

				if (event->mouseButtons() & Qt::LeftButton)
					res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
				else if (event->mouseButtons() & Qt::RightButton)
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