﻿using CodeHelper;
using Controls;
using Controls.Windows;
using ODMR_Lab.IO操作;
using ODMR_Lab.ODMR实验;
using ODMR_Lab.Windows;
using ODMR_Lab.实验类;
using ODMR_Lab.数据处理;
using ODMR_Lab.设备部分;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;
using static System.Windows.Forms.AxHost;
using Label = System.Windows.Controls.Label;

namespace ODMR_Lab
{
    /// <summary>
    /// 实验停止异常
    /// </summary>
    public class ExpStopException : Exception
    {
        public override string Message { get; } = "实验被中止";
    }

    /// <summary>
    /// 实验基类型
    /// </summary>
    /// <typeparam name="ParamType"></typeparam>
    /// <typeparam name="ConfigType"></typeparam>
    public abstract class ExperimentObject<ParamType, ConfigType>
        where ParamType : ExpParamBase
        where ConfigType : ConfigBase
    {

        public ExperimentObject()
        {
            JudgeThreadEndOrResumeAction = JudgeThreadEndOrResume;
        }

        /// <summary>
        /// 非法字符
        /// </summary>
        public static List<char> InvaliidChars = new List<char>() { '$' };

        /// <summary>
        /// 实验文件类型
        /// </summary>
        public abstract ExperimentFileTypes ExpType { get; protected set; }

        /// <summary>
        /// 实验开始时间
        /// </summary>
        public string ExpStartTime { get; set; } = "";

        /// <summary>
        /// 实验结束时间
        /// </summary>
        public string ExpEndTime { get; set; } = "";

        public void SetStartTime(DateTime time)
        {
            ExpStartTime = time.ToString("yyyy-MM-dd HH:mm:ss");
        }

        public void SetEndTime(DateTime time)
        {
            ExpEndTime = time.ToString("yyyy-MM-dd HH:mm:ss");
        }

        /// <summary>
        /// 实验参数
        /// </summary>
        public abstract ParamType Param { get; set; }

        /// <summary>
        /// UI显示参数
        /// </summary>
        public abstract ConfigType Config { get; set; }

        public bool ReadFromFile(string filepath)
        {
            FileObject obj = FileObject.ReadFromFile(filepath);
            if (!obj.Descriptions.Keys.Contains("实验类型"))
            {
                throw new Exception("此文件不属于实验类型文件");
            }

            if (ExpType == (ExperimentFileTypes)Enum.Parse(typeof(ExperimentFileTypes), obj.Descriptions["实验类型"]))
            {
                ExpStartTime = obj.Descriptions["开始时间"];
                ExpEndTime = obj.Descriptions["结束时间"];
                InnerRead(obj);

                if (Param != null)
                {
                    Param.ReadFromDescription(obj.Descriptions);
                }
                if (Config != null)
                {
                    Config.ReadFromDescription(obj.Descriptions);
                }

                return true;
            }
            return false;
        }

        public bool ReadFromExplorer(out string filename)
        {
            filename = "";
            FileObject obj = FileObject.FindFileFromExplorer();
            filename = obj.FilePath;
            if (obj == null) return false;
            if (!obj.Descriptions.Keys.Contains("实验类型"))
            {
                throw new Exception("此文件不属于实验类型文件");
            }

            if (ExpType == (ExperimentFileTypes)Enum.Parse(typeof(ExperimentFileTypes), obj.Descriptions["实验类型"]))
            {
                ExpStartTime = obj.Descriptions["开始时间"];
                ExpEndTime = obj.Descriptions["结束时间"];

                InnerRead(obj);

                if (Param != null)
                {
                    Param.ReadFromDescription(obj.Descriptions);
                }

                if (Config != null)
                {
                    Config.ReadFromDescription(obj.Descriptions);
                }

                return true;
            }
            return false;
        }

        public void WriteToFile(string filefolder, string filename)
        {
            if (!Directory.Exists(filefolder))
            {
                Directory.CreateDirectory(filefolder);
            }

            FileObject obj = new FileObject();

            if (!obj.Descriptions.Keys.Contains("实验类型"))
            {
                obj.Descriptions.Add("实验类型", Enum.GetName(typeof(ExperimentFileTypes), ExpType));
            }
            if (!obj.Descriptions.Keys.Contains("开始时间"))
            {
                obj.Descriptions.Add("开始时间", ExpStartTime);
            }
            if (!obj.Descriptions.Keys.Contains("结束时间"))
            {
                obj.Descriptions.Add("结束时间", ExpEndTime);
            }

            InnerWrite(obj);

            if (Param != null)
            {

                Dictionary<string, string> dic = Param.GenerateDescription();
                foreach (var item in dic)
                {
                    if (!obj.Descriptions.Keys.Contains(item.Key))
                    {
                        obj.Descriptions.Add(item.Key, item.Value);
                    }
                }
            }
            if (Config != null)
            {

                Dictionary<string, string> dic = Config.GenerateDescription();
                foreach (var item in dic)
                {
                    if (!obj.Descriptions.Keys.Contains(item.Key))
                    {
                        obj.Descriptions.Add(item.Key, item.Value);
                    }
                }
            }

            obj.SaveToFile(Path.Combine(filefolder, filename));
        }

        public bool WriteFromExplorer(string defaultName = "")
        {
            FileObject obj = new FileObject();

            if (!obj.Descriptions.Keys.Contains("实验类型"))
            {
                obj.Descriptions.Add("实验类型", Enum.GetName(typeof(ExperimentFileTypes), ExpType));
            }
            if (!obj.Descriptions.Keys.Contains("开始时间"))
            {
                obj.Descriptions.Add("开始时间", ExpStartTime);
            }
            if (!obj.Descriptions.Keys.Contains("结束时间"))
            {
                obj.Descriptions.Add("结束时间", ExpEndTime);
            }

            InnerWrite(obj);

            if (Param != null)
            {
                Dictionary<string, string> dic = Param.GenerateDescription();
                foreach (var item in dic)
                {
                    if (!obj.Descriptions.Keys.Contains(item.Key))
                    {
                        obj.Descriptions.Add(item.Key, item.Value);
                    }
                }
            }

            if (Config != null)
            {
                Dictionary<string, string> dic = Config.GenerateDescription();
                foreach (var item in dic)
                {
                    if (!obj.Descriptions.Keys.Contains(item.Key))
                    {
                        obj.Descriptions.Add(item.Key, item.Value);
                    }
                }
            }

            return obj.SaveFileFromExplorer(defaultname: defaultName);
        }

        /// <summary>
        /// 获取文件的实验类型
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        public static ExperimentFileTypes GetExpType(string filepath)
        {
            Dictionary<string, string> dic = FileObject.ReadDescription(filepath);
            if (!dic.Keys.Contains("实验类型"))
            {
                return ExperimentFileTypes.None;
            }
            try
            {
                return (ExperimentFileTypes)Enum.Parse(typeof(ExperimentFileTypes), dic["实验类型"]);
            }
            catch (Exception ex)
            {
                return ExperimentFileTypes.None;
            }
        }

        /// <summary>
        /// 内部读操作，将obj中的信息转化成对应的ExperimentFileObject
        /// </summary>
        /// <param name="fobj"></param>
        protected abstract void InnerRead(FileObject fobj);

        /// <summary>
        /// 内部写操作，将ExperimentFileObject中的信息转化成FileObject
        /// </summary>
        /// <returns></returns>
        protected abstract void InnerWrite(FileObject obj);

        public DataVisualSource ToDataVisualSource()
        {
            DataVisualSource s = new DataVisualSource();

            s.Params.Add("实验类型", Enum.GetName(ExpType.GetType(), ExpType));
            s.Params.Add("开始时间", ExpStartTime);
            s.Params.Add("结束时间", ExpEndTime);

            InnerToDataVisualSource(s);
            return s;
        }
        /// <summary>
        /// 转换成数据可视化对象
        /// </summary>
        /// <returns></returns>
        protected abstract void InnerToDataVisualSource(DataVisualSource source);

        #region 实验线程部分
        Thread ExpThread { get; set; } = null;

        /// <summary>
        /// 线程运行时的相关控件显示状态
        /// </summary>
        private List<KeyValuePair<FrameworkElement, RunningBehaviours>> ControlStates = new List<KeyValuePair<FrameworkElement, RunningBehaviours>>();

        #region 实验通信控件

        DecoratedButton startButton = null;

        DecoratedButton resumeButton = null;

        DecoratedButton stopButton = null;


        string CurrentexpState = "";
        /// <summary>
        /// 实验进度标签
        /// </summary>
        private TextBlock CurrentexpStateTextBlock = null;

        double CurrentProgress = 0;

        ProgressBar CurrentProgressBar = null;

        private Label ExpStartTimeLabel { get; set; } = null;
        private Label ExpEndTimeLabel { get; set; } = null;

        public void ConnectOuterControl(DecoratedButton StartBtn, DecoratedButton StopBtn, DecoratedButton ResumeBtn, Label StartTimeLabel, Label EndTimeLabel, TextBlock ThreadState, ProgressBar ThreadProgress, List<KeyValuePair<FrameworkElement, RunningBehaviours>> ControlPanels)
        {
            DisConnectOuterControl();
            Dispatcher.CurrentDispatcher.Invoke(() =>
            {
                if (StartBtn != null)
                {
                    startButton = StartBtn;
                    startButton.Click -= StartEvent;
                    startButton.Click += StartEvent;
                }
                if (StopBtn != null)
                {
                    stopButton = StopBtn;
                    stopButton.Click -= StopEvent;
                    stopButton.Click += StopEvent;
                }
                if (ResumeBtn != null)
                {
                    resumeButton = ResumeBtn;
                    resumeButton.Click -= ResumeEvent;
                    resumeButton.Click += ResumeEvent;
                }

                ExpStartTimeLabel = StartTimeLabel;
                if (ExpStartTimeLabel != null)
                {
                    ExpStartTimeLabel.Content = ExpStartTime;
                }
                ExpEndTimeLabel = EndTimeLabel;
                if (ExpEndTimeLabel != null)
                {
                    ExpEndTimeLabel.Content = ExpEndTime;
                }

                CurrentexpStateTextBlock = ThreadState;
                if (CurrentexpStateTextBlock != null)
                {
                    CurrentexpStateTextBlock.Text = CurrentexpState;
                }

                CurrentProgressBar = ThreadProgress;
                if (CurrentProgressBar != null)
                {
                    CurrentProgressBar.Value = CurrentProgress;
                }

                ControlStates = ControlPanels;
                //根据此实验刷新面板状态
                if (IsExpEnd)
                {
                    SetPanelStopState();
                    return;
                }
                if (IsExpResume)
                {
                    SetPanelResumeState();
                    return;
                }
                SetPanelStartState();
            });
        }

        public void DisConnectOuterControl()
        {
            if (startButton != null)
            {
                startButton.Click -= StartEvent;
                startButton = null;
            }
            if (resumeButton != null)
            {
                resumeButton.Click -= ResumeEvent;
                resumeButton = null;
            }
            if (stopButton != null)
            {
                stopButton.Click -= StopEvent;
                stopButton = null;
            }
            if (ExpStartTimeLabel != null)
                ExpStartTimeLabel = null;
            if (ExpEndTimeLabel != null)
                ExpEndTimeLabel = null;
            if (CurrentexpStateTextBlock != null)
                CurrentexpStateTextBlock = null;
            if (CurrentProgressBar != null)
                CurrentProgressBar = null;
            ControlStates = new List<KeyValuePair<FrameworkElement, RunningBehaviours>>();
        }

        public void DisConnectODMRParentExperiment()
        {
            JudgeThreadEndOrResumeAction = JudgeThreadEndOrResume;
        }

        public void ConnectODMRParentExperiment(ODMRExpObject obj)
        {
            JudgeThreadEndOrResumeAction = obj.JudgeThreadEndOrResume;
        }

        /// <summary>
        /// 设置当前线程状态
        /// </summary>
        /// <param name="state"></param>
        public void SetExpState(string state)
        {
            CurrentexpState = state;
            App.Current.Dispatcher.Invoke(() =>
            {
                if (CurrentexpStateTextBlock != null)
                {
                    CurrentexpStateTextBlock.Text = state;
                    CurrentexpStateTextBlock.ToolTip = state;
                }
            });
        }

        /// <summary>
        /// 设置当前线程状态
        /// </summary>
        /// <param name="state"></param>
        public string GetExpState()
        {
            return CurrentexpState;
        }

        /// <summary>
        /// 设置进度(0-100)
        /// </summary>
        public void SetProgress(double value)
        {
            CurrentProgress = value;
            App.Current.Dispatcher.Invoke(() =>
            {
                if (CurrentProgressBar != null)
                {
                    CurrentProgressBar.Value = value;
                }
            });
        }
        #endregion

        /// <summary>
        /// 初始化事件，在读取参数和获取设备后触发
        /// </summary>
        public event Action InitEvent = null;

        public event Action ResumeStateEvent = null;
        public event Action EndStateEvent = null;
        public event Action ErrorStateEvent = null;

        private void SetStartState()
        {
            SetPanelStartState();
        }

        public void SetPanelStartState()
        {
            App.Current.Dispatcher.Invoke(() =>
            {
                TrySetState(startButton, false);
                TrySetState(stopButton, true);
                TrySetState(resumeButton, true);
                foreach (var item in ControlStates)
                {
                    if (item.Value == RunningBehaviours.EnableWhenRunning)
                        item.Key.IsHitTestVisible = true;
                    if (item.Value == RunningBehaviours.DisableWhenRunning || item.Value == RunningBehaviours.DisableWhenRunningEnableWhenResume)
                        item.Key.IsHitTestVisible = false;
                }
            });
            IsExpEnd = false;
            IsExpResume = false;
        }

        private void SetResumeState()
        {
            App.Current.Dispatcher.Invoke(() =>
            {
                SetPanelResumeState();
                ResumeStateEvent?.Invoke();
            });
            IsExpResume = true;
            IsExpEnd = false;
        }

        public void SetPanelResumeState()
        {
            App.Current.Dispatcher.Invoke(() =>
            {
                TrySetState(startButton, true);
                TrySetState(stopButton, true);
                TrySetState(resumeButton, false);
                foreach (var item in ControlStates)
                {
                    if (item.Value == RunningBehaviours.EnableWhenRunning || item.Value == RunningBehaviours.DisableWhenRunningEnableWhenResume)
                        item.Key.IsHitTestVisible = true;
                    if (item.Value == RunningBehaviours.DisableWhenRunning)
                        item.Key.IsHitTestVisible = false;
                }
            });
        }

        private void SetStopState()
        {
            App.Current.Dispatcher.Invoke(() =>
            {
                SetPanelStopState();
                EndStateEvent?.Invoke();
            });
            IsExpEnd = true;
            IsExpResume = false;
        }

        public void SetPanelStopState()
        {
            App.Current.Dispatcher.Invoke(() =>
            {
                TrySetState(startButton, true);
                TrySetState(stopButton, false);
                TrySetState(resumeButton, false);
                foreach (var item in ControlStates)
                {
                    if (item.Value == RunningBehaviours.EnableWhenRunning)
                        item.Key.IsHitTestVisible = false;
                    if (item.Value == RunningBehaviours.DisableWhenRunning || item.Value == RunningBehaviours.DisableWhenRunningEnableWhenResume)
                        item.Key.IsHitTestVisible = true;
                }
            });
        }

        private void TrySetState(FrameworkElement ele, bool state)
        {
            if (ele == null) return;
            else
            {
                if (ele is DecoratedButton)
                {
                    if (state) (ele as DecoratedButton).KeepPressed = false;
                    else
                    {
                        (ele as DecoratedButton).KeepPressed = true;
                    }
                }
                ele.IsEnabled = state;
            }
        }

        public bool IsExpEnd { get; set; } = true;
        public bool IsExpResume { get; set; } = false;

        #region 子实验参数
        /// <summary>
        /// 是否是子实验
        /// </summary>
        public bool IsSubExperiment { get; set; } = false;

        private Exception expFailedException = null;
        public Exception ExpFailedException { get { return expFailedException; } }
        #endregion

        private void StartEvent(object sender, RoutedEventArgs e)
        {
            expFailedException = null;
            if (ThreadResumeFlag == true && ThreadEndFlag == false)
            {
                ThreadResumeFlag = false;
                SetStartState();
                return;
            }
            SetStartState();
            bool IsContinue = true;
            App.Current.Dispatcher.Invoke(() =>
            {
                if (!IsSubExperiment)
                    IsContinue = PreConfirmProcedure();
            });
            if (!IsContinue)
            {
                //设值结束状态
                SetStopState();
                return;
            }
            ThreadEndFlag = false;
            ThreadResumeFlag = false;
            IsExpEnd = false;
            IsExpResume = false;
            ExpThread = new Thread(() =>
            {
                List<InfoBase> Devices = new List<InfoBase>();
                try
                {
                    try
                    {
                        App.Current.Dispatcher.Invoke(() =>
                        {
                            tempConfig = ReadConfig();
                        });
                        //读取参数
                    }
                    catch (Exception ex)
                    {
                        MessageWindow.ShowTipWindow("参数设置存在错误:\n" + ex.Message, MainWindow.Handle);
                        ErrorStateEvent?.Invoke();
                        SetStopState();
                        return;
                    }

                    try
                    {
                        //如果是子实验则不占用设备
                        App.Current.Dispatcher.Invoke(() =>
                        {
                            Devices = GetDevices();
                            DeviceDispatcher.UseDevices(Devices);
                        });
                    }
                    catch (Exception ex)
                    {
                        MessageWindow.ShowTipWindow("设备获取失败:\n" + ex.Message, MainWindow.Handle);
                        ErrorStateEvent?.Invoke();
                        SetStopState();
                        return;
                    }

                    InitEvent?.Invoke();

                    //设置实验时间
                    App.Current.Dispatcher.Invoke(() =>
                    {
                        SetStartTime(DateTime.Now);
                        if (ExpStartTimeLabel != null)
                        {
                            ExpStartTimeLabel.Content = ExpStartTime;
                        }
                        if (ExpEndTimeLabel != null)
                        {
                            ExpEndTimeLabel.Content = "";
                        }
                    });
                    ThreadEndFlag = false;
                    ThreadResumeFlag = false;
                    Config = tempConfig;
                    //进行试验
                    ExperimentEvent();
                    //设值结束状态
                    SetStopState();
                    //结束占用设备
                    DeviceDispatcher.EndUseDevices(Devices);
                    App.Current.Dispatcher.Invoke(() =>
                    {
                        SetEndTime(DateTime.Now);
                        if (ExpEndTimeLabel != null)
                        {
                            ExpEndTimeLabel.Content = ExpEndTime;
                        }
                    });
                }
                catch (Exception ex)
                {
                    DeviceDispatcher.EndUseDevices(Devices);
                    App.Current.Dispatcher.Invoke(() =>
                    {
                        SetEndTime(DateTime.Now);
                        if (ExpEndTimeLabel != null)
                        {
                            ExpEndTimeLabel.Content = ExpEndTime;
                        }
                    });
                    if (ex is ExpStopException)
                    {
                        if (!IsSubExperiment)
                        {
                            MessageWindow.ShowTipWindow("实验已被停止", MainWindow.Handle);
                        }
                    }
                    else
                    {
                        ErrorStateEvent?.Invoke();
                        //如果是子实验
                        if (IsSubExperiment)
                        {
                            expFailedException = ex;
                        }
                        else
                            MessageWindow.ShowTipWindow("实验发生异常,已结束：\n" + ex.Message, MainWindow.Handle);
                    }
                    expFailedException = ex;
                    //设置结束状态
                    SetStopState();
                }
            });
            ExpThread.Start();
        }
        private void ResumeEvent(object sender, RoutedEventArgs e)
        {
            ThreadResumeFlag = true;
            SetExpState("暂停实验...");
        }
        private void StopEvent(object sender, RoutedEventArgs e)
        {
            ThreadEndFlag = true;
            SetExpState("正在停止实验...");
        }

        #region 外部控制
        /// <summary>
        /// 开始实验
        /// </summary>
        public void Start()
        {
            StartEvent(null, new RoutedEventArgs());
        }
        /// <summary>
        /// 暂停实验
        /// </summary>
        public void Resume()
        {
            ResumeEvent(null, new RoutedEventArgs());
        }
        /// <summary>
        /// 停止实验
        /// </summary>
        public void Stop()
        {
            StopEvent(null, new RoutedEventArgs());
        }

        #endregion
        #endregion

        /// <summary>
        /// 线程终止标签
        /// </summary>
        private bool ThreadEndFlag { get; set; } = true;
        /// <summary>
        /// 线程暂停标签
        /// </summary>
        private bool ThreadResumeFlag { get; set; } = false;

        /// <summary>
        /// 线程状态判断函数
        /// </summary>
        public Action JudgeThreadEndOrResumeAction { get; set; } = null;
        /// <summary>
        /// 如果状态为等待则挂起，如果状态为结束则抛出异常
        /// </summary>
        /// <exception cref="Exception"></exception>
        private void JudgeThreadEndOrResume()
        {
            if (ThreadEndFlag)
            {
                throw new ExpStopException();
            }
            if (ThreadResumeFlag)
            {
                SetResumeState();
                while (ThreadResumeFlag)
                {
                    if (ThreadEndFlag)
                    {
                        throw new ExpStopException();
                    }
                    Thread.Sleep(50);
                }
            }
        }



        public void Dispose()
        {
            if (ExpThread == null) return;
            ExpThread.Abort();
            while (ExpThread.ThreadState == ThreadState.Running)
            {
                Thread.Sleep(50);
            }
        }

        /// <summary>
        /// 实验事件
        /// </summary>
        public abstract void ExperimentEvent();


        private ConfigType tempConfig = null;
        /// <summary>
        /// 参数读取事件
        /// </summary>
        public abstract ConfigType ReadConfig();

        /// <summary>
        /// 进行实验前的确认操作,如果不继续则返回false
        /// </summary>
        public abstract bool PreConfirmProcedure();

        /// <summary>
        /// 提供实验需要的设备,返回的对象必须是DeviceInfoBase类
        /// </summary>
        public abstract List<InfoBase> GetDevices();

    }
}
